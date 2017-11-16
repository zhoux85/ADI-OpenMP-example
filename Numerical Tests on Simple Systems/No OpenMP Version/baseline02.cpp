//////////////////////////////////////////////////////////////////////////                                                               
//    Example of FDM Method for 2D heat equation                        
//                                                                      
//          u_t = u_{ xx } +u_{ yy } +f(x, t) 
//                                                                      
//    Test problem : 
//      Exact solution : u(x, y, t) = exp(-t) cos(pi*x) cos(pi*y) 
//      Source term : f(t, x, y) = -u + (2pi^2 - 1)exp(-t)cos(pi*x)cos(pi*y) 
//      (x,y) is in (0,1) and (0,1)                                                         
//     Results :       n              e            ratio               
//                     10           0.0041                              
//     t_final = 0.5   20           0.0010           4.1                
//                     40           2.5192e-04       3.97               
//                     80           6.3069e-05       3.9944             
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
int const nx = 100, ny = 100;//grid numbers
double dx = 1.0 / nx, dy = 1.0 / ny;//space step, 3cm*3cm
double D = 0.001;//D: diffusion coefficient cm^2/ms

/* Time Step */
double dt_max = 1e-05; // Time step (ms)
double dt_min = 1e-06;
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

	for (ncount = 0; ncount <= 8 / dt_max; ncount++){//simulation time is 160ms  160 / dt
		for (i = 1; i <= nx; i++){
			//****no flux boundary conditions*****
			U[i][0] = U[i][1];
			U[i][ny + 1] = U[i][ny];
		}
		for (j = 1; j <= ny; j++){
			U[0][j] = U[1][j];
			U[nx + 1][j] = U[nx][j];
		}

		//*********** step 1 *******
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				Unew[i][j] = U[i][j] + (dt_max / 2)*(U[i + 1][j] + U[i - 1][j] - 2 * U[i][j]) / (dx*dx) \
					+ (dt_max / 2)*(U[i][j + 1] + U[i][j - 1] - 2 * U[i][j]) / (dy*dy);
			}
		}
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				U[i][j] = Unew[i][j];
			}
		}
		//*********** step 1 *******

		//*********** step 2 *******
		dt = dt_max;
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				dUdt[i][j] = usource(U[i][j], i*dx, j*dy, ncount*dt_max);
			}
		}

		int k0, k;
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				// adaptive time step
				if (dUdt[i][j] > 0){
					k0 = 5;
				}else{
					k0 = 1;
				}
				k = k0 + (int)(fabs(dUdt[i][j]) + 0.5); //round the value
				if (k >(int)(dt_max / dt_min)){
					k = (int)(dt_max / dt_min);
				}
				dt = dt_max / k;
				int ttt;
				for (ttt = 0; ttt <= k-1; ttt++){ //from t to t+dt_max, t=t+dt
					U[i][j] = U[i][j] + dt*usource(U[i][j], i*dx, j*dy, ncount*dt_max + ttt*dt);
				}
			}
		}
		//*********** step 2 *******

		//*********** step 3 *******
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				Unew[i][j] = U[i][j] + (dt_max / 2)*(U[i + 1][j] + U[i - 1][j] - 2 * U[i][j]) / (dx*dx) \
					+ (dt_max / 2)*(U[i][j + 1] + U[i][j - 1] - 2 * U[i][j]) / (dy*dy);
			}
		}
		for (i = 1; i <= nx; i++){
			for (j = 1; j <= ny; j++){
				U[i][j] = Unew[i][j];
			}
		}
		//*********** step 3 *******
		t = t + dt_max;
	}

	U[0][0] = U[0][1];
	U[0][ny + 1] = U[0][ny];
	U[nx + 1][0] = U[nx][1];
	U[nx + 1][ny + 1] = U[nx + 1][ny];

	for (i = 0; i <= nx + 1; i++){
		for (j = 0; j <= ny + 1; j++){
			Uexact11[i][j] = uexact(i*dx, j*dy, ncount*dt_max); // Initial Voltage (mv)
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


