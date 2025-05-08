/* 3 state variables -> x, y, vx */
/* controls -> long. acceleration and lateral speed */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>
#include <string.h>
#include "eigen-master/Eigen/Dense"
#include <filesystem>

// include other headers
#include "dionysis_controller.h"
#include "data.h"
#include "controllers.h"
#include "vehicle.h"

using namespace std;

/* Colors */
#define RED   "\x1B[1;31m"
#define GRN   "\x1B[1;32m"
#define RESET "\x1B[0m"

void deleteDirectoryContents(const filesystem::path& dir) {
    for (const auto& entry : filesystem::directory_iterator(dir)) 
        filesystem::remove_all(entry.path());
}
 
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]



// struct Vehicle::Obstacles {
// 	int id;
// 	// double x0, y0, vx0, vy0, ax0, ay0;
// 	vector<double> x, y, vx, vy, ax, ay;
// 	double len, width;
// };

// void Vehicle::IDM(Data data, struct IDM) {
// 	int numObs = data.obs_n, numsteps = data.numsteps;
// 	double T = data.step, ts = data.safety_x; vdx = data.vdx;

// 	double aCon;
//     idm.aMax = 1.5; idm.b = 2.0; idm.d = 4; idm.s0 = 5.0;

// 	init_x.push_back(x0); init_y.push_back(y0); init_vx.push_back(vx0);

// 	for (int k = 0; k < numsteps; k++) {
// 		double dxl = DBL_MAX; 
// 		xl = DBL_MAX; vl = 0.0; hasLeader = 0;
// 		xf = -DBL_MAX; vf = 0.0; hasFollower = 0;
// 		leaderID = -1; followerID = -1;
		
// 		for (int j = 0; j < numObs; j ++){
// 			// TODO: change lane logic

// 			if ((init_y[k] != obs[j].y[k])) continue;
			
// 			dxl = (init_x[k] - obs[j].x[k]);
// 			if (dxl < 0.0 && obs[j].x[k] < xl) {
// 				xl = obs[j].x[k];
// 				vl = obs[j].vx[k];
// 				hasLeader = 1;
// 				leaderID = j;
// 			}

// 			if (dxl > 0.0 && obs[j].x[k] < xl) {
// 				xf = obs[j].x[k];
// 				vf = obs[j].vx[k];
// 				hasFollower = 1;
// 				followerID = j;
// 			}
// 		}
		
// 		/* add MOBIL */
// 		// int laneChange = 0;
// 		// laneChange = mobil(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].leaderID, veh[i].lane[k], veh[i].id, veh[i].followerID, k);
// 		// veh[i].y[k+1] = fmin(veh[i].y[k] + laneChange, noLanes - 1);
// 		init_y.push_back(init_y[k]);

// 		/* calculate longitudinal acceleration */
// 		init_ax.push_back(IDM_accel(data, init_x[k], xl, init_vx[k], vl, vdx, ts, hasLeader));

// 		init_x.push_back(init_x[k] + init_vx[k] * T + 0.5 * init_ax[k] * pow(T, 2));
// 		init_vx.push_back(init_vx[k] + init_ax[k] * T);
		
// 		if (init_vx[k + 1] <= 0.0) {
// 			init_ax[k] = 0.0; 
// 			aCon = idm.b;
// 			init_x[k + 1] = init_x[k] + pow(init_vx[k], 2) / (2.0 * aCon);
// 			init_vx[k + 1] = 0.0;
// 		}
//     }
// }

/*
int Vehicle::mobil(Data data, double x, double xl, double v, double vl, double v0, double ts, int leaderID, double curLane, int vehID, int followerID, int k) {

    int laneChange = 0;
    double Delta_a = 0.5, p = 0.5;
    double delta_a_c = 0.0, delta_a_o = 0.0, delta_a_n = 0.0; 

    if (leaderID == -1) return 0;

    // reference vehicle
    double a_c = IDM_accel(data, x, xl, v, vl, v0, ts, leaderID);
    double a_c_hat = IDM_accel(data, x, xl, v, vl, v0, ts, -1);
    delta_a_c = a_c_hat - a_c;

    // old follower
    double a_o, a_o_hat;
    if (followerID != -1) {
        a_o = IDM_accel(data, veh[followerID].x[k], x, veh[followerID].v[k], v, v0, ts, vehID);
        if (leaderID == -1)
            a_o_hat = IDM_accel(data, veh[followerID].x[k], 0.0, veh[followerID].v[k], 0.0, v0, ts, -1);
        else
            a_o_hat = IDM_accel(data, veh[followerID].x[k], veh[leaderID].x[k], veh[followerID].v[k], veh[leaderID].v[k], v0, ts, leaderID);
        delta_a_o = a_o_hat - a_o;
    }

    find new follower
    int hasNewFollower = 0, newFollowerID;
    for (int i = 0; i < numVeh; i++){
        if (veh[i].id == vehID || veh[i].id == followerID) continue;

        if ( (veh[i].y[k] == veh[vehID].y[k] + 1) && (veh[i].x[k] <= veh[vehID].x[k]) ) {
            hasNewFollower = 1;
            newFollowerID = i;

        }

    }

    if (hasNewFollower == 1) {
        double a_n = idm_accel(x, xl, v, vl, v0, ts, veh[newFollowerID].leaderID);
        double a_n_hat = idm_accel(x, xl, v, vl, v0, ts, vehID);
        delta_a_n = a_n_hat - a_n;

        if (delta_a_n < -3.0) return 0;
    }

    // check if lane change should be applied
    if ( delta_a_c + p*( delta_a_n + delta_a_o ) > Delta_a )
        laneChange = 1;

    return laneChange;
}
*/

// double Vehicle::IDM_accel(Data data, double x, double xl, double v, double vl, double v0, double ts, int leaderID) {

// 	double s = x - xl;
//     double sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
//     double aFree = idm.aMax * (1.0 - pow(v/v0, idm.d));
	
//     double accIDM;
//     if (leaderID == -1)
//         accIDM = aFree;
//     else
//         accIDM = aFree - (idm.aMax * pow(sStar/s, 2));
    
//     return accIDM;

// }

// void Vehicle::initial_path(Data data) {
	
// 	// double uxConst;
// 	// init_x.push_back(x0); init_y.push_back(y0);
// 	// init_vx.push_back(vx0);
// 	// double T = data.step;
	
// 	// uxConst = (data.vdx - data.vx0) / data.numsteps;
// 	// for (int k = 0; k < data.numsteps; k++) {
// 	// 	init_x.push_back(init_x[k] + init_vx[k]*T + 0.5*uxConst*pow(T,2));
// 	// 	init_y.push_back(init_y[k]);
// 	// 	init_vx.push_back(init_vx[k] + uxConst*T);
// 	// }

// 	// TODO: Implement IDM or other method for initial trajectory
// 	// IDM - MOBIL
// 	IDM(data, idm);

// 	fprintf(stderr, "-- Initial path created!\n");
// }

// void Vehicle::obst_prediction(Data data) {

// 	for (int i = 0; i < data.obs_n; i++) {
// 		Obstacles temp_obst;
// 		temp_obst.x.push_back(data.obsx[i]); temp_obst.y.push_back(data.obsy[i]);
// 		temp_obst.vx.push_back(data.obsvx[i]); temp_obst.vy.push_back(data.obsvy[i]);
// 		for (int k = 0; k < data.numsteps; k++) {
// 			temp_obst.x.push_back(temp_obst.x[k] + temp_obst.vx[k]*data.step);
// 			temp_obst.y.push_back(temp_obst.y[k] + temp_obst.vy[k]*data.step);
// 			temp_obst.vx.push_back(temp_obst.vx[k]);
// 			temp_obst.vy.push_back(temp_obst.vy[k]);
// 		}
// 		obs.push_back(temp_obst);
// 	}

// 	fprintf(stderr, "-- Obstacles path created!\n");
// }

int main() {

	clock_t tic, toc;
	double cpu_time;
	
	tic = clock();
	/* start procedure */

	// initialize and read data
	Data input;
	input.read();
	
	// Initialize vehicle(s)
	Vehicle veh(input.x0, input.y0, input.vx0);
	// veh.obst_prediction(input);

	Dionysis_Feedback fbck(1);
	fbck.x0 = input.x0; fbck.y0 = input.y0; fbck.vx0 = input.vx0; fbck.vdx = input.vdx;
	fbck.run();

	fprintf(stderr, "inital_x: %.4f \n", veh.init_x[0]);

return 0;

	// Controller ctr(input.numsteps, input.step, "dddp");
	// ctr.dux0 = 1.0; ctr.duy0 = 1.0;
	// ctr.cx0 = 2.0; ctr.cvx0 = 1.0; ctr.cy0 = 2.0;
	// // ctr.cxm = 5.0 * ctr.dux0; 

	// // ctr.dddp(veh, road, input);
	// ctr.ddp(veh, road, input);

	/* end procedure */
	toc = clock();
	cpu_time = (double)(toc - tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- Total CPU Time used: %.6f --\n", cpu_time);

	return 0;
}
