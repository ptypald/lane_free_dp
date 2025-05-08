#ifndef VEHICLE_H
#define VEHICLE_H

#include <iostream>
#include <vector>

using namespace std;

class Vehicle {
	
	public:
		double x0, y0, vx0;
		double len, width;
		// double tgx, tgy;	// time gaps
		double vdx, vdy;	// desired speeds
		// vehicle bounds
		double UXMIN = -2.0, UXMAX = 1.0;
		double UYMIN = -1.0, UYMAX = 1.0;
		double VXMIN = 0.0, VXMAX = 40.0;

		vector<double> init_vx, init_x, init_y, init_ax;
		vector<double> x, y, vx, ux, uy;

		// IDM parameters
		double xl, vl, xf, vf;
		int leaderID, followerID, hasLeader, hasFollower;

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
		void setDersiredSpeedX(double vx) { vdx = vx; }
		void setDersiredSpeedY(double vy) { vdy = vy; }

	// class structures
	// struct Obstacles;
	// vector<Obstacles> obs;

	// struct IDM{
	// 	double aMax; 
   	// 	double s0;
    // 	double b;
    // 	int d;
	// }idm;

	// class functions
	// void initial_path(Data data);
	// double IDM_accel(Data data, double x, double xl, double v, double vl, double v0, double ts, int leaderID);
	// void IDM(Data data, struct IDM);
	// int mobil(Data data, double x, double xl, double v, double vl, double v0, double ts, int leaderID, double curLane, int vehID, int followerID, int k);
	// void obst_prediction(Data data);
	bool collision_per_veh(double here_x, double obsx, double here_y, double obsy, double here_vx, double obsvx);
};


#endif