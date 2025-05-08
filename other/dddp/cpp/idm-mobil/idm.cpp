#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <vector>

using namespace std; 

#define numVeh 3
#define noLanes 2

double roadWidth = 10.0;

vector<int> obsID { 0, 1 };

struct {
    double T;
    int K;
}P;

struct {
    double aMax; 
    double s0;
    double b;
    int d;
}idm;

struct {
    double *x, *v, *a, *y;
    int id, hasLeader, leaderID;
    int hasFollower, followerID;
    // int *lane;
    double ts, vd, xl, vl, xf, vf;

    // mobil
    int hasNewLeader, hasNewFollower, newLeaderID, newFollowerID;
    double xl_new, vl_new, xf_new, vf_new;
}veh[numVeh];

double idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID) {
    
    double s = x - xl;
    double sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
    double aFree = idm.aMax * (1.0 - pow(v/v0, idm.d));

    double accIDM;
    if (leaderID == -1)
        accIDM = aFree;
    else
        accIDM = aFree - (idm.aMax * pow(sStar/s, 2));
    
    return accIDM;
}

double arrb(double *v, double *a) {

    double vk, ak, Rt, fuel;
    double m = 1600.0, alpha = 0.666, beta1 = 0.0717, beta2 = 0.0344, b1 = 0.269, b2 = 0.0171, b3 = 0.000672;

    fuel = 0.0;
    for (int k = 0; k < P.K+1; k++) {
        vk = v[k]; ak = a[k];
        Rt = b1 + b2 * pow(vk, 2) + b3 * pow(vk, 2) + m * ak / 1000.0;
        if (Rt > 0.0) {
            if (ak > 0.0)
                fuel += alpha + beta1 * Rt * vk + (beta2 * m * pow(ak, 2) * vk / 1000.0);
            else
                fuel += alpha + beta1 * Rt * vk;
        }
        else{
            fuel += alpha;
        }
    }
    // fprintf(stderr, "%.4f \n", fuel * P.T);
    return fuel * P.T;

}

int mobil(double x, double xl, double v, double vl, double v0, double ts, int leaderID, double curLane, int vehID, int followerID, int k) {

    int laneChange = 0;
    double temp_incent = -DBL_MAX;
    double Delta_a = 0.5, p = 0.5;
    double delta_a_c, delta_a_o, delta_a_n;

    if (leaderID == -1) return 0;

    // /* reference vehicle */
    double a_c = idm_accel(x, xl, v, vl, v0, ts, leaderID);

    vector<int> lanes = { 1, -1 };
    for (auto change : lanes) {
        delta_a_c = 0.0, delta_a_o = 0.0, delta_a_n = 0.0; 

        int tragetLane = curLane + change;
        if( tragetLane > noLanes - 1 || tragetLane < 0) continue;;

        // find acceleration in target lane
        double dxl;
        veh[vehID].xl_new = DBL_MAX;
        veh[vehID].hasNewLeader = 0; veh[vehID].newLeaderID = -1;
        for (int ob = 0; ob < numVeh; ob ++){
            // TODO: change lane logic
            if (find(obsID.begin(), obsID.end(), ob) != obsID.end()) continue;
            if (veh[ob].id == veh[vehID].id || veh[ob].id == leaderID) continue;

            if ( (veh[ob].y[k] == curLane + change) && (veh[ob].x[k] > veh[vehID].x[k]) ) {
                dxl = (veh[vehID].x[k] - veh[ob].x[k]);
                if (dxl < 0.0 && veh[ob].x[k] < veh[vehID].xl_new) {
                    veh[vehID].xl_new = veh[ob].x[k];
                    veh[vehID].vl_new = veh[ob].v[k];
                    veh[vehID].hasNewLeader = 1;
                    veh[vehID].newLeaderID = ob;
                }
            }
        }
        double a_c_hat = idm_accel(x, veh[vehID].xl_new, v, veh[vehID].vl_new, v0, ts, veh[vehID].newLeaderID);
        delta_a_c = a_c_hat - a_c;
        
        /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* old follower */
        double a_o, a_o_hat;
        if (followerID != -1) {
            a_o = idm_accel(veh[followerID].x[k], x, veh[followerID].v[k], v, v0, ts, vehID);

            // if the reference vehicle had a leader than the follower adopts him
            if (leaderID == -1)
                a_o_hat = idm_accel(veh[followerID].x[k], 0.0, veh[followerID].v[k], 0.0, v0, ts, -1);
            else
                a_o_hat = idm_accel(veh[followerID].x[k], veh[leaderID].x[k], veh[followerID].v[k], veh[leaderID].v[k], veh[followerID].vd, ts, leaderID);
            delta_a_o = a_o_hat - a_o;
        }

        /* find new follower */
        double dxf;
        veh[vehID].xf_new = -DBL_MAX;
        veh[vehID].hasNewFollower = 0; veh[vehID].newFollowerID = -1;
        for (int j = 0; j < numVeh; j ++){
            // TODO: change lane logic
            if (find(obsID.begin(), obsID.end(), j) != obsID.end()) continue;
            if (veh[j].id == veh[vehID].id || veh[j].id == followerID) continue;

            if ( (veh[j].y[k] == curLane + change) && (veh[j].x[k] > veh[vehID].x[k]) ) {
                dxf = (veh[vehID].x[k] - veh[j].x[k]);
                if (dxf > 0.0 && veh[j].x[k] > veh[vehID].xf_new) {
                    veh[vehID].xf_new = veh[j].x[k];
                    veh[vehID].vf_new = veh[j].v[k];
                    veh[vehID].hasNewFollower = 1;
                    veh[vehID].newFollowerID = j;
                }
            }
        }

        if (veh[vehID].hasNewFollower == 1) {
            double a_n = 0.0;
            double a_n_hat = idm_accel(veh[veh[vehID].newFollowerID].x[k], veh[vehID].x[k], veh[veh[vehID].newFollowerID].v[k], veh[vehID].v[k], veh[veh[vehID].newFollowerID].vd, ts, vehID);
            delta_a_n = a_n_hat - a_n;

            if (delta_a_n < -4.0) return 0;
        }

        /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
        /* check if lane change should be applied */
        double incentive = delta_a_c + p*( delta_a_n + delta_a_o );
        if (vehID == 2)
            fprintf(stderr, "(%d) veh: %d | change: %d | incentive: %.4f \n", k, vehID, change, incentive);
        if ( incentive > Delta_a && incentive > temp_incent) {
            temp_incent = incentive;
            laneChange = change;
        }
    }
    return laneChange;
}

void allocations() {

    for (int n = 0; n < numVeh; n++) {
        veh[n].x = (double*)calloc(P.K + 1, sizeof(double));
        veh[n].y = (double*)calloc(P.K + 1, sizeof(double));
        veh[n].v = (double*)calloc(P.K + 1, sizeof(double));
        veh[n].a = (double*)calloc(P.K + 1, sizeof(double));

        // veh[n].lane = (int*)calloc(P.K + 1, sizeof(int));
    }

}

void exportSimData() {
	
	FILE* fd2 = NULL;
	if ((fd2 = fopen("sim/data/sim.js", "w")) == NULL) {
		printf("\nCouldn't read data file \n");
	}

    int numsteps = P.K;

	double max_pos = 0.0;
	fprintf(fd2, "sim = {\n");
	
    fprintf(fd2, "x: [\n");
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			fprintf(fd2, "%.4f,", veh[idVeh].x[k]);
			if (veh[idVeh].x[k] > max_pos)
				max_pos = veh[idVeh].x[k];
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("y: [\n", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
            double y = (double)veh[idVeh].y[k]*3.3 + 1.5;
            // double y = (double)veh[idVeh].y[k];
			fprintf(fd2, "%.4f,", y);
        }
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vx: [\n", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].v[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("id: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%d,", idVeh + 1);
	fputs("],\n", fd2);

	fputs("Cx: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", numVeh);
	fprintf(fd2, "k:%d,\n", numsteps);
	fprintf(fd2, "roadlength:%f,\n", max_pos);
	fprintf(fd2, "roadwidth:%f,\n", noLanes*3.3);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", P.T);
	fprintf(fd2, "roadbound:%f,\n", 1.2);

	fprintf(fd2, "};\n\n");

	fclose(fd2);
}

void printSol() {
    FILE *fx, *fv, *fa;
    fx = fopen("positions.txt", "w"); fv = fopen("speeds.txt", "w"); fa = fopen("accels.txt", "w");
    for (int i = 0; i < numVeh; i++){
        if (veh[i].id == 0) continue;
        for (int k = 0; k < P.K; k++){
            // fprintf(stderr, "(%d, %d) x: %.4f \t v: %.4f \t a: %.4f \t lane: %d \n", veh[i].id, k, veh[i].x[k], veh[i].v[k], veh[i].a[k], veh[i].lane[k]);
            fprintf(fx, "%.4f \t", veh[i].x[k]); fprintf(fv, "%.4f \t", veh[i].v[k]); fprintf(fa, "%.4f \t", veh[i].a[k]);
        }
        fprintf(fx, "\n"); fprintf(fv, "\n"); fprintf(fa, "\n");
        // double fuel = arrb(veh[i].v, veh[i].a);
        // fprintf(stderr, "fuel: %.4f \n", fuel);
    }

    fclose(fx); fclose(fv); fclose(fa);
}

int main() {
    
    double aCon;
    idm.aMax = 1.5; idm.b = 2.0; idm.d = 4; idm.s0 = 5.0;
    P.T = 0.25; P.K = (int) 30.0/P.T;

    /* allocate necessary memory */
    allocations();

    /* initialize vehicle parameters */
    for (int n = 0; n < numVeh; n++){
        veh[n].ts = 1.0;
        veh[n].id = n;
        veh[n].vd = 16.0;
    }

    /* stopped vehicle static obstacle */
    for (int i : obsID) {
        veh[i].x[0] = 250.0 + 100*i; veh[i].v[0] = 0.0; veh[i].y[0] = (double)i;
    }
    /* initialize vehicles */
    int vehCount = 0;
    for (int i = numVeh - 1; i >= 0; i--) {
        if (find(obsID.begin(), obsID.end(), i) != obsID.end()) continue;
        
        veh[i].x[0] = vehCount * 40.0; veh[i].v[0] = 11.0;
        veh[i].y[0] = 0;
        vehCount++;
    }
    
    for (int i = 0; i < numVeh; i++) {
        fprintf(stderr, "veh: %d | x: %.2f | v: %.2f \n", i, veh[i].x[0], veh[i].v[0]);
    }

    for (int k = 0; k < P.K; k++) {
        for (int i = 0; i < numVeh; i++){
            double dxl = DBL_MAX; 
            veh[i].xl = DBL_MAX; veh[i].vl = 0.0; veh[i].hasLeader = 0;
            veh[i].xf = -DBL_MAX; veh[i].vf = 0.0; veh[i].hasFollower = 0;
            veh[i].leaderID = -1; veh[i].followerID = -1;
            
            if (find(obsID.begin(), obsID.end(), i) != obsID.end()) {
                veh[i].x[k + 1] = veh[i].x[k];
                veh[i].v[k + 1] = veh[i].v[k];
                veh[i].y[k + 1] = veh[i].y[k];
                continue;
            }

            for (int j = 0; j < numVeh; j ++){
                // TODO: change lane logic
                if (veh[j].id == veh[i].id) continue;

                if ((veh[i].y[k] != veh[j].y[k])) continue;

                dxl = (veh[i].x[k] - veh[j].x[k]);
                if (dxl < 0.0 && veh[j].x[k] < veh[i].xl) {
                    veh[i].xl = veh[j].x[k];
                    veh[i].vl = veh[j].v[k];
                    veh[i].hasLeader = 1;
                    veh[i].leaderID = j;
                }

                if (dxl > 0.0 && veh[j].x[k] > veh[i].xf) {
                    veh[i].xf = veh[j].x[k];
                    veh[i].vf = veh[j].v[k];
                    veh[i].hasFollower = 1;
                    veh[i].followerID = j;
                }
            }
            
            /* add MOBIL */
            int laneChange = 0;
            laneChange = mobil(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].leaderID, veh[i].y[k], veh[i].id, veh[i].followerID, k);
            // veh[i].y[k+1] = fmin(veh[i].y[k] + laneChange, noLanes - 1);
            veh[i].y[k+1] = veh[i].y[k] + laneChange;

            /* calculate longitudinal acceleration */
            veh[i].a[k] = idm_accel(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].hasLeader);
            veh[i].x[k + 1] = veh[i].x[k] + veh[i].v[k] * P.T + 0.5 * veh[i].a[k] * pow(P.T, 2);
            veh[i].v[k + 1] = veh[i].v[k] + veh[i].a[k] * P.T;

            if (veh[i].v[k + 1] <= 0.0) {
                veh[i].a[k] = 0.0; aCon = idm.b;
                veh[i].x[k + 1] = veh[i].x[k] + pow(veh[i].v[k], 2) / (2.0 * aCon);
                veh[i].v[k + 1] = 0.0;
            }
        }
    }

    // printSol();
    exportSimData();
    
    return 0;
}