#ifndef DDP_H
#define DDP_H

#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>
// #include "eigen-master/Eigen/Dense"

#include "controllers.h"

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) ((a < b)?(a):(b))
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

using namespace std;

class DDP : public Controllers {
    public:
        int DISPLAY;
    
		vector<double> xmin, xmax, vxmin, vxmax, ymin, ymax;
		vector< vector<double> > X, VX, Y, UXK, UYK;
		vector<double> UX, UY;
		double uxmin, uxmax, uymin, uymax;
        // global and infrastructure specific bounds
        double vxMax, vxMin, yMin, yMax, xMin, xMax;
		
		double opt_J_it, opt_J_prev;
		double ****J, ****Ux, ****Uy;

        DDP() { fprintf(stderr, "-- DDP Initilized with Default Constructor\n"); }
        DDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            fprintf(stderr, "-- DDP Initilized with Parametrized Constructor\n");
        }
        
        void initialize_parameters(Vehicle veh);
        void initial_trajectory();
        void run(Vehicle v);
        int collision(Vehicle v, double x, double y, double vx, int k);
        void allocations();
        void myFree();
        void plotSolution(int it);

        // math functions
        void mat_mult(double **mat1, double **mat2, double **res, int dim1);
        void mat_mult_vect(double **mat1, double *mat2, double *res, int dim1);
        void vect_mult_mat(double *mat1, double **mat2, double *res, int dim1);
        void mat_mult_const(double **mat1, double mat2, double **res, int dim1);
        void vect_mult_const(double *mat1, double mat2, double *res, int dim1);
        void mat_sub(double **mat1, double **mat2, double **res, int dim1);
        void mat_sum(double **mat1, double **mat2, double **res, int dim1);
        void mat_sum2(double *mat1, double *mat2, double *res, int dim1);
        void mat_sub2(double *mat1, double *mat2, double *res, int dim1);
        void mat_2x2_inv(double **mat1, double **res);
        void mat_2x2_transpose(double **mat, double **res);
        double vect_mult_vect(double *mat1, double *mat2, int dim1);
        void vect_mult_vect_to_mat(double *mat1, double *mat2, double **res, int dim1);
};

#endif