#ifndef DDDP_H
#define DDDP_H

#include <iostream>
#include <vector>
#include <math.h>
#include <float.h>

#include "controllers.h"

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) ((a < b)?(a):(b))
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

using namespace std;

class DDDP : public Controllers {
    public:
        int DISPLAY = 0;
        // domain parameters
		vector<int> NX, NVX, NY;
		int maxNX, maxNVX, maxNY;
		vector<double> xmin, xmax, vxmin, vxmax, ymin, ymax;
		vector< vector<double> > X, VX, Y, UXK, UYK;
		vector<double> UX, UY;
		double dx, dvx, dux, dy, duy;
		double uxmin, uxmax, uymin, uymax;
        // global and infrastructure specific bounds
        double vxMax, vxMin, yMin, yMax, xMin, xMax;

		int NUX, NUY;
		// corridor parameters
		double dux0, duy0;
		double cvx0, cxm, cy0, cx0;

		double cx, cy, cvx;
		
		double opt_J_it, opt_J_prev;
		double ****J, ****Ux, ****Uy;

        typedef struct node {
        	int k;
        	int ix, ivx, iy;
        	double x, vx, y;
        	double _x, _vx, _y;
        } node_t;

        typedef struct {
        	double ux, uy;
        } ctrl_t;

        DDDP() { fprintf(stderr, "-- DDDP Initilized with Default Constructor\n"); }
        DDDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            // fprintf(stderr, "-- DDDP Initilized with Parametrized Constructor\n");
        }

        DDDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> uxInit, vector<double> uyInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit, uxInit, uyInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            // fprintf(stderr, "-- DDDP Initilized with Parametrized Constructor\n");
        }

		DDDP(Vehicle veh, int numsteps, double step, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> vyInit, vector<double> uxInit, vector<double> uyInit) : Controllers(veh, numsteps, step, xInit, vxInit, yInit, vyInit, uxInit, uyInit) {
            
            initialize_parameters(veh);
            initial_trajectory();
            
            // fprintf(stderr, "-- DDDP Initilized with Parametrized Constructor\n");
        }
        
        void initialize_parameters(Vehicle veh);
        void initial_trajectory();
        void feasible_domain(int it);
        void run(Vehicle v);
        int collision(Vehicle v, double x, double y, double vx, int k);
        void allocations();
        void myFree();
        int discretizexv(node_t *state, int k);
        void forwardRun(int it);
        void plotDomain(int it);
        void plotSolution(int it);

        node_t *getnext(node_t here, double ux, double uy){
        	static node_t next;
        	int err;

        	next.k = here.k + 1;
        	next.x = next._x = here.x + here.vx*step + 0.5*ux*pow(step, 2);
        	next.y = next._y = here.y + uy*step;
        	next.vx = next._vx = here.vx + ux*step;
            
        	if ((err = discretizexv(&next, next.k))) {
        		fprintf(stderr, " -- INVALID STATE: form (%.3f | %.3f | %.3f) to (%.3f | %.3f | %.3f) at time step %d to %d \n", 
        			here.x, here.y, here.vx, next.x, next.y, next.vx, here.k, next.k);
        		return NULL;
        	}
        	return &next;
        };
};

#endif