#include "idm.h"

void IDM::idm_run() {
	double T = step, ts = 0.2;

	double aCon;
	opt_x.push_back(x0); opt_y.push_back(y0); opt_vx.push_back(vx0); opt_vy.push_back(vy0);
	for (int k = 0; k < numsteps; k++) {
		double dxl = DBL_MAX; 
		xl = DBL_MAX; vl = 0.0; hasLeader = 0;
		xf = -DBL_MAX; vf = 0.0; hasFollower = 0;
		leaderID = -1; followerID = -1;
		
		for (int j = 0; j < numObs; j ++){
			// TODO: change lane logic

			if ((opt_y[k] < obs[j].y[k] + 1.0) &&  (opt_y[k] > obs[j].y[k] - 1.0)) {
			
				dxl = (opt_x[k] - obs[j].x[k]);
				if (dxl < 0.0 && obs[j].x[k] < xl) {
					xl = obs[j].x[k];
					vl = obs[j].vx[k];
					hasLeader = 1;
					leaderID = j;
				}

				if (dxl > 0.0 && obs[j].x[k] < xl) {
					xf = obs[j].x[k];
					vf = obs[j].vx[k];
					hasFollower = 1;
					followerID = j;
				}
			}
		}
		
		/* add MOBIL */
		// int laneChange = 0;
		// laneChange = mobil(veh[i].x[k], veh[i].xl, veh[i].v[k], veh[i].vl, veh[i].vd, veh[i].ts, veh[i].leaderID, veh[i].lane[k], veh[i].id, veh[i].followerID, k);
		// veh[i].y[k+1] = fmin(veh[i].y[k] + laneChange, noLanes - 1);
		opt_y.push_back(opt_y[k]);

		/* calculate longitudinal acceleration */
		opt_ux.push_back(idm_accel(opt_x[k], xl, opt_vx[k], vl, vdx, ts, hasLeader));
		opt_uy.push_back(0.0);

		opt_x.push_back(opt_x[k] + opt_vx[k] * T + 0.5 * opt_ux[k] * pow(T, 2));
		opt_vx.push_back(opt_vx[k] + opt_ux[k] * T);
		opt_vy.push_back(0.0);
		
		if (opt_vx[k + 1] <= 0.0) {
			opt_ux[k] = 0.0; 
			aCon = b;
			opt_x[k + 1] = opt_x[k] + pow(opt_vx[k], 2) / (2.0 * aCon);
			opt_vx[k + 1] = 0.0;
		}
    }
	plotSolution(opt_x, opt_y, opt_vx, opt_ux);
}

double IDM::idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID) {

	double s = x - xl;
    double sStar = s0 + fmax(v*ts + ((v * (v - vl)) / (2.0 * sqrt(aMax * b))), 0.0);
    double aFree = aMax * (1.0 - pow(v/v0, d));
	
    double accIDM;
    if (leaderID == -1)
        accIDM = aFree;
    else
        accIDM = aFree - (aMax * pow(sStar/s, 2));
    
    return accIDM;

}

void IDM::plotSolution(vector<double> opt_x, vector<double> opt_y, vector<double> opt_vx, vector<double> opt_ux) {
	char fileBuffer[64];
    snprintf(fileBuffer, sizeof(char) * 64, "../outputs/idm/solution/solution_%02i.py", 0);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");

	/* python file */
	fputs("solution = {\n", fout);
	
	// vx
	fputs("'vx': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_vx[k]);
	fputs("],\n", fout);

	fputs("'vy': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", 0.0);
	fputs("],\n", fout);

	// x
	fputs("'x': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_x[k]);
	fputs("],\n", fout);

	// y
	fputs("'y': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_y[k]);
	fputs("],\n", fout);

	// ux
	fputs("'ux': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_ux[k]);
	fputs("],\n", fout);

	fputs("'uy': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", 0.0);
	fputs("],\n", fout);

	fprintf(fout, "'k':%d,\n", numsteps);
	fprintf(fout, "'Step':%f,\n", step);
	fprintf(fout, "}\n\n");


	fclose(fout);
};