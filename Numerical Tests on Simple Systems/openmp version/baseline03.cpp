//////////////////////////////////////////////////////////////////////////                                                               
//    Example of FDM Method for 2D heat equation                        
//                                                                      
//          u_t = u_{ xx } +u_{ yy } +f(x, t) 
//                                                                      
//    Test problem : 
//      Exact solution : u(x, y, t) = exp(-t) cos(pi*x) cos(pi*y) 
//      Source term : f(t, x, y) = -u + (2pi^2 - 1)exp(-t)cos(pi*x)cos(pi*y) 
//      (x,y) is in (0,1) and (0,1)                                                           
//////////////////////////////////////////////////////////////////////////
/* modified by Xiang Zhou, 2017/11/5 */

//#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>  
#include <sys/timeb.h>

#define PI 3.141592

//*******FDM parameters for LR91 *******
int const nx = 10, ny = 10;//grid numbers
double dx = 1.0 / nx, dy = 1.0 / ny;//space step, 3cm*3cm

/* Time Step */
double dt_max = (1e-05); // Time step (ms)
double dt_min = (1e-06);
double dt; // Time step (ms)
double t; // Time (ms)
int steps; // Number of Steps
int increment; // Loop Control Variable

/* Voltage */
double U[nx + 2][ny + 2]; // Initial Voltage (mv)
double Unew[nx + 2][ny + 2]; // Initial Voltage (mv)
double Uerror[nx + 2][ny + 2]; //  error 

FILE *single_ap;
void performance();
double uexact(double x, double y, double t);
double usource(double uu, double x, double y, double t);

int main(int argc, char* argv[])
{
	/* Data File */
	FILE *fevaluation;
	fevaluation = fopen("fevaluation", "w");

	/* Time Loop Conditions */
	t = 0.0; // Time (ms)

	int ncount, i, j;

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			U[i][j] = uexact(i*dx, j*dy, 0); // Initial Voltage (mv)
		}
	}

	struct timeb start, end;
	int diff;
	ftime(&start);

	double Uexact11[nx + 2][ny + 2] = { 0 };
	double dUdt[nx + 2][ny + 2] = { 0 };
	steps = 8 / dt_max;
	omp_set_num_threads(8);
	int tid;

	// for ADI of step 1 and 3
	double belta[nx + 1];
	double eta = dt_max / (dx*dx);
	double b = 1 + eta;
	double b_1 = 1 + eta/2;//take care the boundary value
	double b_n = 1 + eta/2;//take care the boundary value
	double c = -eta / 2;
	double a = -eta / 2;
	double f[nx + 1][ny + 1];
	double y_temp[nx + 1];

#pragma omp parallel private(ncount,tid)
	{
		tid = omp_get_thread_num();// Get the thread ID
		printf("ddd from thread %d\n", tid);
		for (ncount = 1; ncount <= steps; ncount++){//simulation time is 160ms  160 / dt
			if (tid == 0){
				for (i = 1; i < nx + 1; i++){
					//****no flux boundary conditions*****
					U[i][0] = U[i][1];
					U[i][ny + 1] = U[i][ny];
				}
				for (j = 1; j < ny + 1; j++){
					U[0][j] = U[1][j];
					U[nx + 1][j] = U[nx][j];
				}
			}
#pragma omp barrier

			//*********** step 1 *******
#pragma omp for private(i,j) schedule(static)
			for (j = 1; j < ny + 1; j++){
				for (i = 1; i < nx + 1; i++){
					if (j == 1){
						f[i][j] = U[i][j] + (eta / 2)*(U[i][j] - 2 * U[i][j] + U[i][j + 1]);
					}else if (j == ny){
						f[i][j] = U[i][j] + (eta / 2)*(U[i][j - 1] - 2 * U[i][j] + U[i][j]);
					}else{
						f[i][j] = U[i][j] + (eta / 2)*(U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]);
					}
				}
			}
#pragma omp barrier

#pragma omp for private(i,j,belta,y_temp) schedule(static)
			for (j = 1; j < ny + 1; j++){
				belta[1] = c / b_1;
				y_temp[1] = f[1][j] / b_1;
				for (i = 2; i < nx; i++){ //i = 2,3,...,n-1
					belta[i] = c / (b - a*belta[i - 1]);
					y_temp[i] = (f[i][j] - a*y_temp[i - 1]) / (b - a*belta[i - 1]);
				}
				y_temp[nx] = (f[nx][j] - a*y_temp[nx - 1]) / (b_n - a*belta[nx - 1]);
				U[nx][j] = y_temp[nx];
				for (i = nx - 1; i >= 1; i--){
					U[i][j] = y_temp[i] - belta[i] * U[i + 1][j];
				}
			}
#pragma omp barrier
			//*********** step 1 *******
			//if (tid==0){
			//	printf("aaaa");
			//}
			//*********** step 2 *******
			dt = dt_max;
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i <= nx; i++){
				for (j = 1; j <= ny; j++){
					dUdt[i][j] = usource(U[i][j], i*dx, j*dy, (ncount - 1)*dt_max);
				}
			}
#pragma omp barrier

			int k0, k;
#pragma omp for private(i,j, k0, k, dt) schedule(static)
			for (i = 1; i <= nx; i++){
				for (j = 1; j <= ny; j++){
					// adaptive time step
					if (dUdt[i][j] > 0){
						k0 = 5;
					}else{
						k0 = 1;
					}
					k = k0 + (int)(fabs(dUdt[i][j]) + 0.5); //round the value
					if (k > (int)(dt_max / dt_min)){
						k = (int)(dt_max / dt_min);
					}
					dt = dt_max / k;
					int ttt;
					for (ttt = 0; ttt <= k - 1; ttt++){ //from t to t+dt_max, t=t+dt
						U[i][j] = U[i][j] + dt*usource(U[i][j], i*dx, j*dy, (ncount - 1)*dt_max + ttt*dt);
					}
				}
			}
#pragma omp barrier
			//*********** step 2 *******
			//*********** step 3 *******
#pragma omp for private(i,j) schedule(static)
			for (i = 1; i < nx + 1; i++){
				for (j = 1; j < ny + 1; j++){
					if (i == 1){
						f[i][j] = U[i][j] + (eta / 2)*(U[i][j] - 2 * U[i][j] + U[i + 1][j]);
					}else if (i == nx){
						f[i][j] = U[i][j] + (eta / 2)*(U[i - 1][j] - 2 * U[i][j] + U[i][j]);
					}else{
						f[i][j] = U[i][j] + (eta / 2)*(U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]);
					}
				}
			}
#pragma omp barrier

#pragma omp for private(i,j,belta,y_temp) schedule(static)
			for (i = 1; i < nx + 1; i++){
				belta[1] = c / b_1;
				y_temp[1] = f[i][1] / b_1;
				for (j = 2; j < ny; j++){
					belta[j] = c / (b - a*belta[j - 1]);
					y_temp[j] = (f[i][j] - a*y_temp[j - 1]) / (b - a*belta[j - 1]);
				}
				y_temp[ny] = (f[i][ny] - a*y_temp[ny - 1]) / (b_n - a*belta[ny - 1]);
				U[i][ny] = y_temp[ny];
				for (j = ny - 1; j >= 1; j--){
					U[i][j] = y_temp[j] - belta[j] * U[i][j + 1];
				}
			}
#pragma omp barrier
			//*********** step 3 *******
			t = t + dt_max;
		}
	}

	U[0][0] = U[0][1];
	U[0][ny + 1] = U[0][ny];
	U[nx + 1][0] = U[nx][1];
	U[nx + 1][ny + 1] = U[nx + 1][ny];

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			Uexact11[i][j] = uexact(i*dx, j*dy, steps*dt_max); // Initial Voltage (mv)
		}
	}

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			Uerror[i][j] = fabs(Uexact11[i][j] - U[i][j]); // Initial Voltage (mv)
		}
	}

	//matrix norm L_1 error
	double sum_err = 0, err_L1 = 0;
	for (j = 0; j <= ny + 1; j++){
		for (i = 0; i <= nx + 1; i++){
			sum_err += Uerror[i][j];
		}
	}
	err_L1 = sum_err*dx*dy;

	ftime(&end);
	diff = (int)(1000.0*(end.time - start.time) + (end.millitm - start.millitm));
	fprintf(fevaluation, "%d s\n%g\n", diff / 1000, err_L1);
	fclose(fevaluation);
}

double uexact(double x, double y, double t){
	return  exp(-t)*cos(PI*x)*cos(PI*y);
}

double usource(double uu, double x, double y, double t){
	return -uu + (2 * PI*PI)*exp(-t)*cos(PI*x)*cos(PI*y);
}


