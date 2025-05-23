#ifndef DDP_H
#define DDP_H

#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>
#include "eigen-master/Eigen/Dense"
 
#include "controllers.h"

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) ((a < b)?(a):(b))
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

using namespace std;

class DDP : public Controllers {
    public:
        int DISPLAY;

        // boundary controller parameters
        double Klat = 1.0;
    
		vector<double> xmin, xmax, vxmin, vxmax, ymin, ymax;
		vector< vector<double> > X, VX, Y, UXK, UYK;
		vector<double> UX, UY;
		double uxmin, uxmax, uymin, uymax;
        // global and infrastructure specific bounds
        double vxMax, vxMin, yMin, yMax, xMin, xMax;
		
		double opt_J_it, opt_J_prev;
		double ****J, ****Ux, ****Uy;

        // ellipsoid parameters
        double p = 5.0, gamma = 2.0;

        DDP() { fprintf(stderr, "-- DDP Initilized with Default Constructor\n"); }
        
        DDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            fprintf(stderr, "-- DDP Initilized with Parametrized Constructor\n");
        }

        DDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> uxInit, vector<double> uyInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit, uxInit, uyInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            fprintf(stderr, "-- DDP Initilized with Parametrized Constructor\n");
        }

        DDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> vyInit, vector<double> uxInit, vector<double> uyInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit, vyInit, uxInit, uyInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            fprintf(stderr, "-- DDP Initilized with Parametrized Constructor\n");
        }

        typedef struct Stages {
            Eigen::MatrixXf alphaStore;
	        Eigen::MatrixXf betaStore;
            Eigen::MatrixXf PStore;
	        Eigen::MatrixXf QStore;
        }stages;
        
        void initialize_parameters(Vehicle veh);
        void initial_trajectory();
        void run(Vehicle v);
        // void plotSolution(int it);
        void plotSolution(vector<double> opt_x, vector<double> opt_y, vector<double> opt_vx, vector<double> opt_vy, vector<double> opt_ux, vector<double> opt_uy, int it);

        double axUp(double axPoint); double axLow(double axPoint);
        double ayUp(double ayPoint); double ayLow(double ayPoint);
        double ayLB(double yPoint, double vyPoint, double ayPoint, double dy, double dvy, double day);
        double ayUB(double yPoint, double vyPoint, double ayPoint, double dy, double dvy, double day);
        double testUB(double dy, double dvy, double yPoint, double vyPoint, double ayPoint);
        double testLB(double dy, double dvy, double yPoint, double vyPoint, double ayPoint);
        vector<int> activeSet(vector<int> S0, double xP, double yP, double vxP, double vyP, double axP, double ayP, int k);
        vector<int> modifyS(vector<int> S0, Eigen::MatrixXf p, Eigen::MatrixXf lambda, double axP, double ayP, double vyP, double yP, int k);
};

#endif