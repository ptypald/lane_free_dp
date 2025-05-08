#include "vehicle.h"

void Vehicle::obst_prediction() {

	for (int i = 0; i < obs_n; i++) {
		Obstacles temp_obst;
		temp_obst.x.push_back(obsx0[i]); temp_obst.y.push_back(obsy0[i]);
		temp_obst.vx.push_back(obsvx0[i]); temp_obst.vy.push_back(obsvy0[i]);
		for (int k = 0; k < numsteps; k++) {
			temp_obst.x.push_back(temp_obst.x[k] + temp_obst.vx[k]*step);
			temp_obst.y.push_back(temp_obst.y[k] + temp_obst.vy[k]*step);
			temp_obst.vx.push_back(temp_obst.vx[k]);
			temp_obst.vy.push_back(temp_obst.vy[k]);
		}
		obs.push_back(temp_obst);
	}

	fprintf(stderr, "-- Obstacles path created!\n");
}