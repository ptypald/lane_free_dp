#ifndef DATA_H
#define DATA_H
#include <iostream>
#include <vector>
#include <string.h>
#include <vector>

using namespace std;

class Data {
	public:
		// flags
		int DISPLAY;
		// input
		int numsteps;
		double step, timeHorizon;
		double safety_x, safety_y;
		double x0, y0, vx0, ax0, ay0; // is this needed?
		double vdx, vdy;
		vector<double> obsx, obsy, obsvx, obsvy, obsax, obsay;
		int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len;
		int obs_n;

		// road data
		double xMin, xMax, yMin, yMax;

		// output
		vector<double> x, y, vx, ux, uy; // is this needed?

		// functions
		void read();
};

#endif