#ifndef VEHICLE_H
#define VEHICLE_H

#include <vector>
#include "data.h"
using namespace std;

class Vehicle {
	
	private:
		double x0, y0, vx0;
		double len, width;
		double vdx, vdy;	// desired speeds
		// vehicle bounds
		double UXMIN, UXMAX, UYMIN, UYMAX, VXMIN, VXMAX;

	public:
		vector<double> init_vx, init_x, init_y;
		vector<double> x, y, vx, ux, uy;
		double safety_x, safety_y;	// time gaps
		int numsteps;
		double step;

		int obs_n;
		vector<double> obsx0, obsy0, obsvx0, obsvy0;
		struct Obstacles {
			int id;
			vector<double> x, y, vx, vy, ax, ay;
			double len, width;
		};
		vector<Obstacles> obs;

		// Parametrized Constructor
		Vehicle(double x, double y, double vx) { 
			x0 = x; y0 = y; vx0 = vx;
		}
		
		// getter functions
		double getX0() { return x0; }				double getY0() { return y0; }
		double getVX0() { return vx0; }
		double getLength() { return len; } 			double getWidth() { return width; }
		double getDersiredSpeedX() { return vdx; } 	double getDersiredSpeedY() { return vdy; }
		double getUXmin() { return UXMIN; } 		double getUXmax() { return UXMAX; }
		double getUYmin() { return UYMIN; } 		double getUYmax() { return UYMAX; }
		double getVXmin() { return VXMIN; } 		double getVXmax() { return VXMAX; }
		//setter functions
		void setDersiredSpeedX(double vx) { vdx = vx; }		void setDersiredSpeedY(double vy) { vdy = vy; }
		void setUXmin(double uxmin) { UXMIN = uxmin; }		void setUXmax(double uxmax) { UXMAX = uxmax; }
		void setUYmin(double uymin) { UYMIN = uymin; }		void setUYmax(double uymax) { UYMAX = uymax; }
		void setVXmin(double vxmin) { VXMIN = vxmin; }		void setVXmax(double vxmax) { VXMAX = vxmax; }

		void obst_prediction();
};		

#endif