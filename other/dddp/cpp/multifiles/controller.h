#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <vector>
#include <iostream>
#include "vehicle.h"

using namespace std;

class Controller {

	public:
		// infrastructure parameters
		double roadxMin, roadxMax, roadyMin, roadyMax;
		// vehicle capabilities
		double uxUB, uxLB, uyUB, uyLB, vxUB, vxLB;
		// general parameters
		int numsteps;
		int max_iter;
		double step;
		string method;
		// domain parameters
		vector<int> NX, NVX, NY;
		int maxNX, maxNVX, maxNY;
		vector<double> xmin, xmax, vxmin, vxmax, ymin, ymax;
		vector<vector<double> > X, VX, Y;
		vector<double> UX, UY;
		double dx, dvx, dux, dy, duy;
		double uxmin, uxmax, uymin, uymax;
		int NUX, NUY;
		// vehicle info
		double x0, y0, vx0, vdx, vdy;
		
		// optimal trajectories
        vector<double> opt_x, opt_vx, opt_y;
        double opt_J_it, opt_J_prev;
		double ****J, ****Ux, ****Uy;		

	// parametrized constructor
	Controller(int K, double T, string m) {
		numsteps = K;
		step = T;
		method = m;
	}
	Controller() { }

	typedef struct node {
		int k;
		int ix, ivx, iy;
		double x, vx, y;
		double _x, _vx, _y;
	} node_t;

	typedef struct {
		double ux, uy;
	} ctrl_t;

	void myAlloc();
	void myFree();
	int discretizexv(node_t *state, int k);
	int collision(Vehicle v, double here_x, double here_y, double here_vx, int here_k);
	void printSolution(Vehicle v, int it);
	// void plotDomain(Vehicle v, int it);

	node_t *getnext(node_t here, double ux, double uy) {
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