#ifndef D_CONTROLLER_H
// #define EXTERN_C
// #else
// #define EXTERN_C extern
// #endif
#define D_CONTROLLER_H

#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

/* Colors */
#define RED   "\x1B[1;31m"
#define GRN   "\x1B[1;32m"
#define RESET "\x1B[0m"

class Dionysis_Feedback {
	public:
        double x0, y0, vx0, vdx;
        vector<double> obsx0, obsy0, obsvx0;

		/* Parameters */
		double a = 5.0;           	// perpendicular of the road in meters (half road width)
		double sigma = 5.0;			// Vehicle length in meters
		double phi = 0.25;        	// Angle max
		// #define p 5.11251       	// 1/(3*pow(tan(phi),2))
		double p = 15.0;       		// 1/(3*pow(tan(phi),2))
		double L = 5.5;       		// sigma*fmax(2*sqrt(p)*sin(phi),sqrt(1+(p-1)*pow(sin(phi),2)))

		int N;						// time step

		double epsilon = 0.2;
		double lambda = 25;       	// radius of communication
		double c = 1.5;          	// 3.86866 //1.5
		double vmax = 35;         	// Maximum velocity of vehicles

		int A = 1;
		int n;             			// number of vehicles(+1)

		int pot = 0;
		int model;

		double m1, m2, qu, q;

		// VISOSCITY      
		//double visc = 0.05/pow(vmax,2);      
		// double visc =  1 / vmax;
		double visc = 0.0;   		// 0 if we dont want to take other vehicle speeds into account      

		int tmax = 20; 				// Time ending point         
		double T = 0.00125; 			// Step-size  
		int flag_fail = 0;

		int every = 1000;

		int counT = 1;
		double Fmax = 0;

	Dionysis_Feedback(int m) {
		model = m;
		if (m == 1) {
			// MODEL NEWTONIAN
			m1 = 0.5;
			// #define m2   1/vmax
			m2 = 0.5;
			qu = 1;
			q = 0.003;
		}
		else if (m == 2) {
			// MODEL RELATIVISTIC
			m1 = 0.3;
			m2 = 1/pow(vmax,2);
			qu =  1;
			q = 0.000003;
		}
		else {
			fprintf(stderr, RED "Error:" RESET "unknown model \n");
			exit(1);
		}
		
	}

	double distance(double xi, double xj, double yi, double yj, double var);
	double Vfunc(double var1, int var2, double var3, double var4);
	double sum(double(*), int var2);
	double func1(double var2, double var1);
	double Ufunc(double var1, double var2, double var3, int var4);
	double kappa(double var1, double var2, double var3, double var4);
	int mod(int v1, int v2);
	void run();
};


#endif