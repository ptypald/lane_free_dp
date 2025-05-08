#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdbool.h>

#define numVeh 2
#define noLanes 2

struct {
    double T;
    int K;
    int simulationTime;
}P;

struct {
    double aMax; 
    double s0;
    double b;
    int d;
}idm;

struct {
    double *x, *v, *a;
    int id, hasLeader, *leaderID;
    int hasFollower, followerID;
    int *lane;
    double ts, vd, xl, vl, xf, vf;
}veh[numVeh];

double idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID) {
    
    double s = x - xl;
    double sStar = idm.s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(idm.aMax * idm.b))), 0.0);
    double aFree = idm.aMax * (1.0 - pow(v/v0, idm.d));

    double accIDM;
    if (leaderID == -1)
        accIDM = aFree;
    else if (leaderID >= 0)
        accIDM = aFree - (idm.aMax * pow(sStar/s, 2));
    else if (leaderID == -2){
        // accIDM = aFree;
        // accIDM = idm.aMax * (1.0 - pow(x/220.0, 4));
        // accIDM = idm.aMax * ( (1.0 - pow(x/220.0, 4)) + (1.0 - pow(v/v0, 6)) ) / 2.0;
        accIDM = idm.aMax * fmin( (1.0 - pow(x/220.0, 4)), (1.0 - pow(v/v0, 2)) );
    }
    else
        fprintf(stderr, "Unknown Leader Type. \n");
    
    
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
    double Delta_a = 0.5, p = 0.5;
    double delta_a_c = 0.0, delta_a_o = 0.0, delta_a_n = 0.0; 

    if (leaderID == -1 || leaderID == 0) return 0;

    // /* reference vehicle */
    double a_c = idm_accel(x, xl, v, vl, v0, ts, leaderID);
    double a_c_hat = idm_accel(x, xl, v, vl, v0, ts, -1);
    delta_a_c = a_c_hat - a_c;

    // /* old follower */
    double a_o, a_o_hat;
    if (followerID != -1) {
        a_o = idm_accel(veh[followerID].x[k], x, veh[followerID].v[k], v, v0, ts, vehID);
        if (leaderID == -1)
            a_o_hat = idm_accel(veh[followerID].x[k], 0.0, veh[followerID].v[k], 0.0, v0, ts, -1);
        else
            a_o_hat = idm_accel(veh[followerID].x[k], veh[leaderID].x[k], veh[followerID].v[k], veh[leaderID].v[k], v0, ts, leaderID);
        delta_a_o = a_o_hat - a_o;
    }

    /* DO LOOP FOR EACH LANE */
    for (int i = 0; i < noLanes; i++) {
        
    }
    // /* new follower */
    int hasNewFollower = 0;
    if (hasNewFollower == 1) {
        double a_n = idm_accel(x, xl, v, vl, v0, ts, vehID);
        double a_n_hat = idm_accel(x, xl, v, vl, v0, ts, leaderID);
        delta_a_n = a_n_hat - a_n;
    }

    /* check if lane change should be applied */
    if ( delta_a_c + p*( delta_a_n + delta_a_o ) > Delta_a )
        laneChange = 1;

    return laneChange;
}

void allocations() {

    for (int n = 0; n < numVeh; n++) {
        veh[n].x = (double*)calloc(1024 + 1, sizeof(double));
        veh[n].v = (double*)calloc(1024 + 1, sizeof(double));
        veh[n].a = (double*)calloc(1024 + 1, sizeof(double));

        veh[n].lane = (int*)calloc(1024 + 1, sizeof(int));
        veh[n].leaderID = (int*)calloc(1024 + 1, sizeof(int));
    }

}

int main() {
    
    double aCon;
    idm.aMax = 3.0; idm.b = 3.0; idm.d = 4; idm.s0 = 5.0;
    P.T = 0.25; P.K = (int) 30.0/P.T;

    /* allocate necessary memory */
    allocations();

    /* initialize vehicle parameters */
    for (int n = 0; n < numVeh; n++){
        veh[n].ts = 1.0;
        veh[n].id = n;
        veh[n].vd = 16.0;
    }

    /* stopped vehicle representing traffic signal */
    veh[0].x[0] = 150.0; veh[0].v[0] = 0.0;
    /* initialize vehicles */
    int vehCount = 0;
    for (int i = 0; i < numVeh; i++) {
        if (i == 0) continue;
        veh[i].x[0] = vehCount * 50.0; veh[i].v[0] = 11.0;
        vehCount++;
    }
    
    int k = 0;
    while(true) {
        
        // if (k >= P.K - 1) {
        if (veh[1].x[k] >= 220.0) {
            P.simulationTime = k;
            break;
        }
        for (int i = 0; i < numVeh; i++){
            if (veh[i].id == 0) {
                veh[i].x[k + 1] = veh[i].x[k];
                veh[i].v[k + 1] = veh[i].v[k];
                continue;
            }

            double dxl = DBL_MAX; 
            veh[i].xl = DBL_MAX; veh[i].vl = 0.0; veh[i].hasLeader = 0;
            veh[i].xf = -DBL_MAX; veh[i].vf = 0.0; veh[i].hasFollower = 0;
            veh[i].leaderID[k] = -1; veh[i].followerID = -1;
            for (int j = 0; j < numVeh; j ++){
                if (veh[j].id == veh[i].id || (veh[j].lane[k] != veh[i].lane[k] && j != 0)) continue;

                dxl = (veh[i].x[k] - veh[j].x[k]);
                if (dxl < 0.0 && veh[j].x[k] < veh[i].xl) {
                    veh[i].xl = veh[j].x[k];
                    veh[i].vl = veh[j].v[k];
                    veh[i].hasLeader = 1;
                    veh[i].leaderID[k] = j;
                }

                if (dxl > 0.0 && veh[j].x[k] < veh[i].xl) {
                    veh[i].xf = veh[j].x[k];
                    veh[i].vf = veh[j].v[k];
                    veh[i].hasFollower = 1;
                    veh[i].followerID = j;
                }
            }
            
            // if (k >= (int)(30/P.T) && veh[i].leaderID[k] == 0) {
            if (k >= (int)(15/P.T)) {
                veh[i].leaderID[k] = -2;
                veh[i].vd = 11.0;
            }

            /* add MOBIL */
            int laneChange = mobil(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].leaderID[k], veh[i].lane[k], veh[i].id, veh[i].followerID, k);
            veh[i].lane[k+1] = fmin(veh[i].lane[k] + laneChange, noLanes - 1);

            /* calculate longitudinal acceleration */
            veh[i].a[k] = idm_accel(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].leaderID[k]);
            veh[i].x[k + 1] = veh[i].x[k] + veh[i].v[k] * P.T + 0.5 * veh[i].a[k] * pow(P.T, 2);
            veh[i].v[k + 1] = veh[i].v[k] + veh[i].a[k] * P.T;

            if (veh[i].v[k + 1] <= 0.0) {
                veh[i].a[k] = 0.0; aCon = idm.b;
                veh[i].x[k + 1] = veh[i].x[k] + pow(veh[i].v[k], 2) / (2.0 * aCon);
                veh[i].v[k + 1] = 0.0;
            }
        }

        k++;
    }

    FILE *fx, *fv, *fa;
    fx = fopen("positions.txt", "w"); fv = fopen("speeds.txt", "w"); fa = fopen("accels.txt", "w");
    for (int i = 0; i < numVeh; i++){
        if (veh[i].id == 0) continue;
        for (int k = 0; k <= P.simulationTime; k++){
            fprintf(stderr, "(%d, %d) x: %.4f \t v: %.4f \t a: %.4f \t lane: %d \t leader: %d \n", veh[i].id, k, veh[i].x[k], veh[i].v[k], veh[i].a[k], veh[i].lane[k], veh[i].leaderID[k]);
            fprintf(fx, "%.4f \t", veh[i].x[k]); fprintf(fv, "%.4f \t", veh[i].v[k]); fprintf(fa, "%.4f \t", veh[i].a[k]);
        }
        fprintf(fx, "\n"); fprintf(fv, "\n"); fprintf(fa, "\n");
        double fuel = arrb(veh[i].v, veh[i].a);
        fprintf(stderr, "fuel: %.4f \n", fuel);
    }
    
    return 0;
}