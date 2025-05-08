#ifndef DATA_H
#define DATA_H

#include <vector>

using namespace std;
 
class Data {
	public:
		// flags
		int DISPLAY;
		// input
		int numsteps;
		double step;
		double safety_x, safety_y;
		double x0, y0, vx0, ax0, ay0; // is this needed?
        double UXMIN, UXMAX, UYMIN, UYMAX, VXMIN, VXMAX;
		double vdx, vdy;
		vector<double> obsx, obsy, obsvx, obsvy, obsax, obsay;
		int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len;
		int obs_n;

		// output
		vector<double> x, y, vx, ux, uy; // is this needed?

		// functions
		void input();
		void exportSimData(class Vehicle veh, class Infrastructure inf, int it);
		// void exportPlotData(Vehicle veh);
};

#endif