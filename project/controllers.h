#ifndef CONTROLLERS_H
#define CONTROLLERS_H

#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>

#include "vehicle.h"

using namespace std;

class Controllers {

	public:
		// general parameters
		int numsteps;
		double step;
		string method;

		// vehicles parameters
		int numVeh;		// number of controlled and obstacle vehicles respectively
		double tgx, tgy;		// long. and lateral time-gaps
		double vdx;

		// obstacles parameters
		int numObs;
		// vector<double> obsx0, obsvx0, obsy0, obsvy0;

		// optimal trajectories
		vector<double> opt_x, opt_vx, opt_y, opt_vy, opt_ux, opt_uy;
		// initial trajectories
		double x0, y0, vx0, vy0, ax0, ay0;
		vector<double> init_x, init_vx, init_y, init_vy, init_ux, init_uy;
		
		// class structures
		struct Obstacles {
			int id;
			vector<double> x, y, vx, vy, ax, ay;
			double len, width;
		};
		vector<Obstacles> obs;

		// parametrized constructor
		Controllers(){ fprintf(stderr, "-- Controllers Initilized with Default Constructor\n"); };
		Controllers(Vehicle veh, int K, double T) 
		: numsteps(K), step(T), vdx(veh.vdx), numObs(veh.numObs),
		  x0(veh.x0), y0(veh.y0), vx0(veh.vx0), obs(numObs) {

			obst_prediction(veh.obsx0, veh.obsy0, veh.obsvx0);
			fprintf(stderr, "-- Controllers Initilized with Parametrized Constructor\n");
		}

		Controllers(Vehicle veh, int K, double T, vector<double> xInit, vector<double> vxInit, vector<double> yInit) 
		: numsteps(K), step(T), vdx(veh.vdx), numObs(veh.numObs),
		  x0(veh.x0), y0(veh.y0), vx0(veh.vx0), 
		  init_x(xInit), init_vx(vxInit), init_y(yInit),
		  obs(numObs) {

			obst_prediction(veh.obsx0, veh.obsy0, veh.obsvx0);
			fprintf(stderr, "-- Controllers Initilized with Parametrized Constructor with Initial Trajectory\n");
		}

		Controllers(Vehicle veh, int K, double T, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> uxInit, vector<double> uyInit) 
		: numsteps(K), step(T), vdx(veh.vdx), numObs(veh.numObs),
		  x0(veh.x0), y0(veh.y0), vx0(veh.vx0), 
		  init_x(xInit), init_vx(vxInit), init_y(yInit), init_ux(uxInit), init_uy(uyInit),
		  obs(numObs) {

			obst_prediction(veh.obsx0, veh.obsy0, veh.obsvx0);
			fprintf(stderr, "-- Controllers Initilized with Parametrized Constructor with Initial Trajectory\n");
		}

		Controllers(Vehicle veh, int K, double T, vector<double> xInit, vector<double> vxInit, vector<double> yInit, vector<double> vyInit, vector<double> uxInit, vector<double> uyInit) 
		: numsteps(K), step(T), vdx(veh.vdx), numObs(veh.numObs),
		  x0(veh.x0), y0(veh.y0), vx0(veh.vx0), 
		  init_x(xInit), init_vx(vxInit), init_y(yInit), init_vy(vyInit), init_ux(uxInit), init_uy(uyInit),
		  obs(numObs) {

			obst_prediction(veh.obsx0, veh.obsy0, veh.obsvx0);
			fprintf(stderr, "-- Controllers Initilized with Parametrized Constructor with Initial Trajectory\n");
		}

		// Class Functions
		void obst_prediction(vector<double> x0, vector<double> y0, vector<double> vx0);
		void exportSimData(int it);
};


#endif