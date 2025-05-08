#ifndef DATA_H
#define DATA_H

#include <vector>

using namespace std;

class Data {
	public:
		// flags
		int DISPLAY = 0;
		// input
		int numsteps = 24; // 6 seconds
		double step = 1.0;
		double safety_x, safety_y;
		double x0, y0, vx0, ax0, ay0; // is this needed?
		double vdx, vdy;
		vector<double> obsx, obsy, obsvx, obsvy, obsax, obsay;
		int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len;
		int obs_n;

		// road data
		double xMin = 0.0, xMax = 100000.0, yMin = 0.0, yMax = 10.0;

		// output
		vector<double> x, y, vx, ux, uy; // is this needed?

		// functions
		void read();
		// void exportSimData(class Vehicle veh, class Infrastructure inf, int it, int numsteps);
		// void exportPlotData(Vehicle veh);
};

#endif