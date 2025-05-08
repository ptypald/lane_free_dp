#ifndef VEHICLE_H
#define VEHICLE_H

#include <iostream>
#include <vector>

#include "data.h"

using namespace std;

class Vehicle {
	
	public:
		double x0, y0, vx0;
		double len, width;
		int numObs;
		// double tgx, tgy;	// time gaps
		double vdx, vdy;	// desired speeds
		// vehicle bounds
		double UXMIN, UXMAX, UYMIN, UYMAX, VXMIN, VXMAX;

		vector<double> obsx0, obsy0, obsvx0, obsvy0, obsax0, obsay0;
		vector<double> init_vx, init_x, init_y, init_ax;
		vector<double> x, y, vx, ux, uy;

		// Parametrized Constructor
		Vehicle() : UXMIN(-2.0), UXMAX(1.0), UYMIN(-1.0), UYMAX(1.0), VXMIN(0.0), VXMAX(40.0) {};
		Vehicle(Data input) : UXMIN(-2.0), UXMAX(1.0), UYMIN(-1.0), UYMAX(1.0), VXMIN(0.0), VXMAX(40.0) 
		{ 
			x0 = input.x0; y0 = input.y0; vx0 = input.vx0;
			vdx = input.vdx; numObs = input.obs_n;
			for (int i = 0; i < numObs; i++) {
				obsx0.push_back(input.obsx[i]); 
				obsy0.push_back(input.obsy[i]);
				obsvx0.push_back(input.obsvx[i]);
			}
		}

		// class functions
		int collision_per_obs(double here_x, double obsx, double here_y, double obsy, double here_vx, double obsvx);
};


#endif