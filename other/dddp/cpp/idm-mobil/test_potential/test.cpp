#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

// #define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
// #include <crtdbg.h>
#define sscanf_s sscanf

#define MAX(a,b) (a > b)?(a):(b)
#define MIN(a,b) (a < b)?(a):(b)
#define NU 2
#define NX 5

#define Power(a, b) pow(a, b)
#define Sqrt(a) sqrt(a)

#define SMOOTHMAX(a, b) 0.5 * (a + b + sqrt(pow(a - b, 2) + EPSILON))
#define SMOOTHMIN(a, b) SMOOTHMAX(-a, -b)
#define DSM(a, b) 0.5 * (-1.0 - ((-a + b)/(sqrt(EPSILON + pow((-a + b), 2)))))

#define CONSTRAINTS 1
#define SDIR 3

#define DP 0
#define EPSILON 1.0

double r2 = 1.3;
double l_m = 1.4;
double cmu = 1.5;
double smu = 0.95;

/* ellipse parameters */
double p1 = 24, p2 = 24;
double eTerm = 1.0;

// 	  ux,  uy,  vd,  off, ellipse_soft, vy
enum {				   PUX, PUY,  PVD,  POF,  PCR,  PAX, PNV, PVY };
double penalty[] = { 1.0, 1.0, 3.0, 15.0, 125.0, 1.0, 1.0, 1.0 }; //trajectories
double c_min[] = { -12.0, -2.0 };	// { -24, -4 }
double c_max[] = { 8.0, 2.0 };		// { -16,  4 }
#define AXMIN -6
#define AXMAX 4

double **states, **costates, **x_opt;
double **control, **u_opt;
double **reduced_g;
double **r, r1, **shift;

double  **sum_term, **sum_term_2;
double **states_prev, **states_cur, **control_prev, **control_cur, **costates_cur, **reduced_g_cur;
double **states_a1, **costates_a1, **reduced_g_a1;
double **states_a2, **costates_a2, **reduced_g_a2;
double **control_a1, **control_a2;
double **states_alpha, **control_alpha, **costates_alpha, **reduced_g_alpha;
double **grad_dif;
double **reduced_g_prev, **s_dir, **temp_states, **temp_costates, **temp_control, **temp_reduced_g, **control_pre;

const char *criterion = "gradient";
struct{
	int iter = 1500;		/* max number of iterations */
	double accur = 0.5;	/* terminal criterion accuracy */
	/* FDS parameters */
	double a = 0.1;		/* line search initial step*/
	double sigma = 0.1;
	double rho = 0.01;
	double tau = 1.2;
	int lin_it = 100;
	int restart_it = 10;
}AP;

struct {
	double vd, T, C;
	double x0, y0, v0, a0;
	double *xk, *yk, *vk;
	double *ax_init;

	double **obst_x, **obst_y, **obst_v;
	int **obst_connected;
	double **obst_ax;
	double *obst_planTime, *obst_vd;
	int *obst_id;
	int obst_x_len, obst_y_len, obst_v_len;

	double *ux0, *uy0;
	int ux0_len, uy0_len;

	int numlanes, numsteps, n;

	double safety;
	double planTime;
	int id;
	int init_crash;
}P = { 0 };

struct {
	double **u;
	double **x;
}init = {0};

struct {
	double ***x;
	double **v;

	/* for connected */
	double **ax;
	/*double **y;*/
}obs = { 0 };

static double inner_product(double **A, double **B, int numsteps, int nu) {
	int i, j;
	double result;
	result = 0.0;

	for (i = 0; i < nu; i++){
		for (j = 0; j < numsteps; j++){
			result += A[i][j] * B[i][j];
		}
	}
	return result;
}

static double sign(double a){
	if (a > 0.0){
		return 1.0;
	}
	else if (a < 0.0){
		return -1.0;
	}else{
		return 0.0;
	}
}

static bool restart(double **s_dir, double **reduced_g, int numsteps, int nu){

	if(inner_product(s_dir, reduced_g, numsteps, NU) >= -0.001*sqrt(inner_product(s_dir, s_dir, numsteps, NU)*inner_product(reduced_g, reduced_g, numsteps, NU))){
		
		return true;
	}
	return false;
}

static double quadratic_interpolation(double a, double b, double dfa, double dfb){
	double alpha;

	alpha = (b*dfa - a*dfb)/(dfa-dfb);

	return alpha;
}

static double cubic_interpolation(double a, double b, double fa, double fb, double dfa, double dfb) {
	double z, w;
	double alpha;

	z = 3.0* ((fa - fb) / (b - a)) + dfa + dfb;
	w = pow(z, 2) - dfa * dfb;

	if (w > 0) {
		w = sqrt(w);
		alpha = b - (((dfb + w - z) / (dfb - dfa + 2.0*w)) * (b - a));
	}
	else {
		fprintf(stderr, "negative sqrt in cubic interpolation, apply bisection \n");
		alpha = (a + b) / 2.0;
	}
	return alpha;
}

static bool termination(const char *termination, double accuracy, double cost_cur, double cost_prev, double **gradient, int numsteps){
	double temp;

	
	if (strncmp (termination, "relative", strlen(termination)) == 0){
		temp = (cost_cur - cost_prev) / (cost_cur + 0.001);
		if (fabs(temp) < accuracy){
			return true;
		}
	}
	else if (strncmp (termination, "gradient", strlen(termination)) == 0){
		temp = inner_product(gradient, gradient, numsteps, NU);
		if (temp < accuracy){
			return true;
		}
	}
	return false;
}

static void adjustConnectedObsPaths(int numsteps, double T, int obsInd, int off) {
	int j;
	double *temp_x, *temp_y, *temp_v;
	assert((temp_x = (double*)calloc(sizeof(double), numsteps + 1)));
	assert((temp_y = (double*)calloc(sizeof(double), numsteps + 1)));
	assert((temp_v = (double*)calloc(sizeof(double), numsteps + 1)));

	if (off <= P.numsteps - 1) {
		for (j = 0; j < numsteps - off; j++) {
			obs.x[obsInd][0][j + 1] = obs.x[obsInd][0][j] + T * obs.v[obsInd][j] + pow(T, 2) * P.obst_ax[obsInd][j + off];
			obs.x[obsInd][1][j + 1] = ((3.0 * P.numlanes) - 1.5) - (3.0 * P.obst_y[obsInd][j + 1 + off]);
			obs.v[obsInd][j + 1] = obs.v[obsInd][j] + T * P.obst_ax[obsInd][j + off];
		}

		for (j = numsteps - off; j < numsteps; j++) {
			obs.x[obsInd][0][j + 1] = obs.x[obsInd][0][j] + T * obs.v[obsInd][numsteps - off];
			obs.x[obsInd][1][j] = obs.x[obsInd][1][numsteps - off - 1];
			obs.v[obsInd][j + 1] = obs.v[obsInd][numsteps - off];
		}
	}
	else {
		temp_x[0] = P.obst_x[obsInd][0]; temp_y[0] = P.obst_y[obsInd][0]; temp_v[0] = P.obst_v[obsInd][0];
		for (j = 0; j < numsteps; j++) {
			temp_x[j + 1] = temp_x[j] + T * temp_v[j] + pow(T, 2) * P.obst_ax[obsInd][j];
			temp_y[j + 1] = ((3.0 * P.numlanes) - 1.5) - (3.0 * P.obst_y[obsInd][j + 1]);
			temp_v[j + 1] = temp_v[j] + T * P.obst_ax[obsInd][j];
			//fprintf(stderr, "tempx: %.3f -- tempy: %.3f -- tempv: %.3f \n", temp_x[j], temp_y[j], temp_v[j]);
		}

		obs.x[obsInd][0][j] = temp_x[numsteps - 1]; //temp_x[numsteps-1];
		for (j = 0; j < numsteps; j++) {
			obs.x[obsInd][0][j + 1] = obs.x[obsInd][0][j] + T * temp_v[numsteps-1];
			obs.x[obsInd][1][j] = temp_y[numsteps-1];
			obs.v[obsInd][j + 1] = temp_v[numsteps-1];
		}

	}


}

/* new */
static void connectedInfeasiblePaths(int numsteps, double T, int obsInd) {
	int j;

	double egoX = P.x0;
	double egoV = P.v0;
	int temp = 0;

	for (j = 0; j < numsteps; j++) {
		egoX = egoX + egoV * T;	// +0.5*P.a0*pow(T, 2);
		//egoV = egoV + P.a0*T;
		double discreteEgo = P.y0;
		double discreteObs = round((((3.0 * P.numlanes) - 1.5) - (obs.x[obsInd][1][j])) / 3.0);
		//fprintf(stderr, "(%d) dx %.3f -- obsty: %.2f (%.2f)\n", P.obst_id[obsInd], abs(obs.x[obsInd][0][j] - egoX), discreteObs, discreteEgo);
		//if ((abs(obs.x[obsInd][0][j] - egoX) < 5.5*P.C) && (discreteEgo == discreteObs)) {
		if ( ( (abs(obs.x[obsInd][0][j] - egoX) < (P.v0*P.safety + P.C) / (2.0)) ) && (discreteEgo == discreteObs)) {
			temp = 1;
			fprintf(stderr, "Infeasible lane change: %d \n", P.obst_id[obsInd]);
			break;
		}
	}
	if (temp == 1) {
		/* if Infeasible Lnae Change then the obstacle turns to zombie */
		obs.x[obsInd][0][0] = P.obst_x[obsInd][0];
		obs.x[obsInd][1][0] = ((3.0 * P.numlanes) - 1.5) - (3.0 * P.obst_y[obsInd][0]);
		obs.v[obsInd][0] = P.obst_v[obsInd][0];
		for (j = 0; j < numsteps; j++) {
			obs.x[obsInd][0][j + 1] = obs.x[obsInd][0][j] + T * obs.v[obsInd][0];
			obs.x[obsInd][1][j + 1] = obs.x[obsInd][1][j];
			obs.v[obsInd][j + 1] = obs.v[obsInd][0];
		}
		//fprintf(stderr, "obst %d  turned to zombie\n", P.obst_id[obsInd]);
	}

}

static void obst_trajectories(int numsteps, double T) {
	int i, j;
	double L = cmu * P.C;

	assert((obs.x = (double***)calloc(sizeof(double**), P.n)));
	assert((obs.v = (double**)calloc(sizeof(double*), P.n)));
	assert((obs.ax = (double**)calloc(sizeof(double*), P.n)));
	assert((r = (double**)calloc(sizeof(double*), P.n)));
	assert((shift = (double**)calloc(sizeof(double*), P.n)));
	for (i = 0; i < P.n; i++) {
		assert((obs.x[i] = (double**)calloc(sizeof(double*), NU)));
		assert((obs.v[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((obs.ax[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((r[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((shift[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		for (int j = 0; j < NU; j++)
			assert((obs.x[i][j] = (double*)calloc(sizeof(double), numsteps + 1)));
	}	

	for (i = 0; i < P.n; i++) {
		obs.x[i][0][0] = P.obst_x[i][0];
		obs.x[i][1][0] = ((3.0 * P.numlanes) - 1.5) - (3.0 * P.obst_y[i][0]);
		//obs.x[i][1][0] = ((3.0 * P.obst_y[i][0]) + 1.5);
		obs.v[i][0] = P.obst_v[i][0];

		int off = (P.obst_planTime[i] == 0.0) ? 0 : int((P.planTime - P.obst_planTime[i]) / P.T);
		//fprintf(stderr, "obs: %d -- off: %d \n", P.obst_id[i], off);
		if (off < 0)
			fprintf(stderr, "Error: Negative 'off' (%lf)\n", off);
		if ((P.obst_connected[i][0] == 1.0) && ( off == 0 || P.obst_planTime[i] == 0)) {
		/* connected vehicles with full trajectory */
			for (j = 0; j < numsteps; j++) {
				obs.x[i][0][j + 1] = obs.x[i][0][j] + T * obs.v[i][j] + pow(T, 2) * P.obst_ax[i][j];
				obs.x[i][1][j + 1] = ((3.0 * P.numlanes) - 1.5) - (3.0 * P.obst_y[i][j + 1]);
				obs.v[i][j + 1] = obs.v[i][j] + T * P.obst_ax[i][j];
			}
			connectedInfeasiblePaths(numsteps, T, i);
			//continue;
		}
		else if ((P.obst_connected[i][0] == 1.0) && (off > 0 && P.obst_planTime[i] > 0)) {
			/* connected vehicles which trajectories need to be cropped in time */
			//fprintf(stderr, "adjust the path of obs %d with off: %d \n", P.obst_id[i], off);
			adjustConnectedObsPaths(numsteps, T, i, off);
		
			/* new - if there is obstacle infeasible path, turn the connected obstacle to zombie */
			connectedInfeasiblePaths(numsteps, T, i);
		}
		else {
		/* no connected vehicles (zombies) */
			for (j = 0; j < numsteps; j++) {
				obs.x[i][0][j + 1] = obs.x[i][0][j] + T * obs.v[i][0];
				obs.x[i][1][j + 1] = obs.x[i][1][j];
				obs.v[i][j + 1] = obs.v[i][0];
			}
		}
	}
}

static void compute_states(double **x_init, double **control, int numsteps, double T, double **states_calc){
	int i;

	for(i = 0; i < NX; i++){
		states_calc[i][0] = x_init[i][0];
	}

	for (i = 0; i < numsteps; i++){
			states_calc[0][i + 1] = states_calc[0][i] + states_calc[2][i] * T + 0.5*states_calc[4][i] * pow(T, 2) + (1.0 / 6.0)*control[0][i] * pow(T, 3);
			states_calc[1][i + 1] = states_calc[1][i] + states_calc[3][i] * T + 0.5 * control[1][i] * pow(T, 2);
			states_calc[2][i + 1] = states_calc[2][i] + states_calc[4][i] * T + 0.5*control[0][i] * pow(T, 2);
			states_calc[3][i + 1] = states_calc[3][i] + control[1][i] * T;
			states_calc[4][i + 1] = states_calc[4][i] + control[0][i] * T;
	}
}

static void compute_costates(double **states, double **control, int numsteps, double T, double **costates_calc){
	int i, j;
	double x1new, x2new, off_road;
	double dvdv, dvn, dax, dvy;
	double y0 = (P.numlanes * 3.0) - l_m;
	double L = cmu * P.C;

	for (i = 0; i < NX; i++) {
		costates_calc[i][numsteps] = 0.0;	/* all costates equal to zero because we dont have terminal conditions */
	}

	i = numsteps - 1;
	while (i>=0){

		double costate_term_1 = 0.0;
		double costate_term_2 = 0.0;
		double costate_term_3 = 0.0;

		/* collision avoidance term derivatice w.r.t x, y, v */
		for (j = 0; j < P.n; j++) {
			// x1new = obs.x[j][0][i] - states[0][i];
			// x2new = obs.x[j][1][i] - states[1][i];

			double obsv = obs.v[j][i], obsx = obs.x[j][0][i], obsy = obs.x[j][1][i];
			double x = states[0][i], y = states[1][i], v = states[2][i];
			// double safety = P.safety;
			double dx = x - obsx;
			double dy = y - obsy;

			// costate_term_1 += -((Power(2, p1)*eTerm*p1*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), -1 + p1)) /
			// 	((2 * L + obsv * safety + safety * v)*Power(1 + Power(2, p1)*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), p1) + Power((obsy - y) / r2, p2), 2)));
			// costate_term_2 += (eTerm*p2*Power((obsy - y) / r2, -1 + p2)) / (r2*Power(1 + Power(2, p1)*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), p1) + Power((obsy - y) / r2, p2), 2));
			// costate_term_3 += -((Power(2, p1)*eTerm*p1*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), -1 + p1)*(safety / (2.*(2 * L + obsv * safety + safety * v)) - (safety*(-obsx + (-(obsv*safety) + safety * v) / 2. + x)) / Power(2 * L + obsv * safety + safety * v, 2))) /
			// 	Power(1 + Power(2, p1)*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), p1) + Power((obsy - y) / r2, p2), 2));

			costate_term_1 += 0.5*((-0.05*(0. + x))/Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)) - (0.05*(0. + x)*(1 - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2))))/(Sqrt(0.1 + Power(1 - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)),2))*Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2))));
			costate_term_2 += 0.5*((-0.5*(-5. + y))/Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)) - (0.5*(1 - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)))*(-5. + y))/(Sqrt(0.1 + Power(1 - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)),2))*Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2))));

		}

		/* road boundaries term derivative w.r.t y */
		if (states[1][i] < l_m) { 
			off_road = penalty[3] * ((states[1][i] - l_m) / pow(y0, 2));
		}
		else if (states[1][i] > y0) {
			off_road = penalty[3] * ((states[1][i] - y0) / pow(y0, 2));
		}
		else
			off_road = 0.0;

		dvdv = penalty[2] * ((states[2][i] - P.vd) / pow(P.vd, 2));
		dvy = penalty[7] * ((states[3][i]) / (c_max[1])) * (1.0 / (c_max[1]));
		dax = penalty[5] * ((states[4][i]) / (AXMAX)) * (1.0 / (AXMAX));
		dvn = penalty[6] * (SMOOTHMIN(states[2][i], 0.0) / P.vd) * DSM(states[2][i], 0.0) / P.vd;

		/* costates calculation */
		costates_calc[0][i] = costates_calc[0][i + 1] + penalty[4] * costate_term_1;
		costates_calc[1][i] = costates_calc[1][i + 1] + penalty[4] * costate_term_2 + off_road;
		costates_calc[2][i] = costates_calc[2][i + 1] + T * costates_calc[0][i + 1] + dvdv + penalty[4] * costate_term_3 + dvn;
		costates_calc[3][i] = costates_calc[3][i + 1] + T * costates_calc[1][i + 1] + dvy; //new
		costates_calc[4][i] = costates_calc[4][i + 1] + T * costates_calc[2][i + 1] + 0.5 * pow(T, 2) * costates[0][i + 1] + dax;
		
		i = i - 1;
	}
}

static double compute_cost(double **states, double **control, int numsteps, double T){
	int i, j;
	double ux, uy, vd, vn, ax, elipse_cost, vy;
	double dx, dy, off_road_cost;
	double y0 = (P.numlanes * 3.0) - l_m;
	double cost = 0.0;
	double L = cmu * P.C;
	double rEgo, rObs;

	for (i = 0; i < numsteps; i++){
		elipse_cost = 0.0;

		/* collision avoidance term */
		for (j = 0; j < P.n; j++) {
			double obsv = obs.v[j][i], obsx = obs.x[j][0][i], obsy = obs.x[j][1][i], x = states[0][i], y = states[1][i], v = states[2][i], safety = P.safety;

			dx = x - obsx;
			dy = y - obsy;
			
			// rEgo = (v * safety) + L;
			// rObs = (obsv * safety) + L;
			// r[j][i] = (rEgo + rObs) / 2.0;
			// shift[j][i] = (rObs - rEgo) / 2.0;

			// elipse_cost += eTerm / (1 + Power(2, p1)*Power((-obsx + (-(obsv*safety) + safety * v) / 2. + x) / (2 * L + obsv * safety + safety * v), p1) + Power((dy) / r2, p2));
			elipse_cost += 0.5*(1 + Sqrt(0.1 + Power(1 - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)),2)) - Sqrt(0.05*Power(0. + x,2) + 0.5*Power(-5. + y,2)));
		}
		elipse_cost *= penalty[4];

		/* terms for u (longtitudinal/lateral) and desired speed */
		double norm_ux = (control[0][i]) / (c_min[0]), norm_uy = (control[1][i]) / (c_min[1]);
		double norm_vd = (states[2][i] - P.vd) / P.vd;
		double norm_vy = (states[3][i]) / (c_max[1]);
		double norm_ax = (states[4][i]) / (AXMAX);
		double norm_vn = SMOOTHMIN(states[2][i], 0.0) / P.vd;

		ux = 0.5 * penalty[0] * pow(norm_ux, 2);
		uy = 0.5 * penalty[1] * pow(norm_uy, 2);
		vd = 0.5 * penalty[2] * pow(norm_vd, 2);
		vy = 0.5 * penalty[7] * pow(norm_vy, 2);
		ax = 0.5 * penalty[5] * pow(norm_ax, 2);
		vn = 0.5 * penalty[6] * pow(norm_vn, 2);

		/* road boundaries term */
		if (states[1][i] < l_m) {
			off_road_cost = 0.5 * penalty[3] * pow((states[1][i] - l_m) / y0, 2);
		}
		else if (states[1][i] > y0) {
			off_road_cost = 0.5 * penalty[3] * pow((states[1][i] - y0) / y0, 2);
		}
		else {
			off_road_cost = 0.0;
		}

		cost += ux + uy + vd + vy + off_road_cost + elipse_cost + ax + vn;
		
	}

	return cost;
}

static void compute_reduced_gradient(double **states, double **costates, double **control, int numsteps, double T, double **reduced_g_calc){
	int i, j;

	for (i = 0; i < numsteps; i++) {
		double dnorm_ux = ((control[0][i]) / (c_min[0])) * (1.0 / (c_min[0]));
		double dnorm_uy = ((control[1][i]) / (c_min[1])) * (1.0 / (c_min[1]));

		reduced_g_calc[0][i] = penalty[0] * dnorm_ux + 0.5*pow(T, 2)*costates[2][i + 1] + T * costates[4][i + 1] + (1.0 / 6.0)*pow(T, 3)*costates[0][i + 1];
		reduced_g_calc[1][i] = penalty[1] * dnorm_uy + 0.5*pow(T, 2)*costates[1][i + 1] + T * costates[3][i + 1];
	}
	if (CONSTRAINTS){
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				if( (control[i][j] == c_min[i]) && (reduced_g_calc[i][j] > 0.0) ){
					reduced_g_calc[i][j] = 0.0;
				}
				else if( (control[i][j] == c_max[i]) && (reduced_g_calc[i][j] < 0.0) ){
					reduced_g_calc[i][j] = 0.0;
				}
			}
		}
	}
}

static void control_bound(double **control, double *control_min, double *control_max, int numsteps, int nu){
	int i, j;

	for (i = 0; i < nu; i++){
		for (j = 0; j < numsteps; j++){
			if(control[i][j] < control_min[i]){
				control[i][j] = control_min[i];
			}
			else if(control[i][j] > control_max[i]){
				control[i][j] = control_max[i];
			}
		}
	}
}

static double line_search(double a, double **states, double **control, double **s_dir, double sigma, double tau, double rho, int lin_it, int numsteps, double T){
	int i, j, b_it;
	double f0, df0, fai_1, fai, dfai;
	double a_cur, a_prev, a_temp, a1 = 0.0, a2 = 0.0;
	double fa1, dfa1, fa2, dfa2;
	double fa, dfa, alpha;

	f0 = compute_cost(states, control, numsteps, T);
	compute_costates(states, control, numsteps, T, costates);
	compute_reduced_gradient(states, costates, control, numsteps, T, reduced_g);
	df0 = inner_product(reduced_g, s_dir, numsteps, NU);

	// if (fabs(df0) < 10e-8){
	// 	printf("Line search exits - no possibility dF/d_alpha = %f\n", df0);
	// 	return 0;
	// }

	a_cur  = a;		/* a_cur: initial step, given from the user */
	a_prev = 0.0;	/* Initialization of a_prev */

	/* Bracketing */
	/*            */
	for (b_it = 0; b_it < lin_it; b_it++){
		/* Calculation of previous control*/
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				control_prev[i][j] = control[i][j] + a_prev*s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(control_prev, c_min, c_max, numsteps, NU);
		}
		/* Calculation of previous states */
		compute_states(states, control_prev, numsteps, T, states_prev);
		/* Calculation of F(a_i-1) */
		fai_1 = compute_cost(states_prev, control_prev, numsteps, T);

		/* Calculation of current control */
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				control_cur[i][j] = control[i][j] + a_cur*s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(control_cur, c_min, c_max, numsteps, NU);
		}

		/* Calculation of current states */
		compute_states(states, control_cur, numsteps, T, states_cur);
		/* Calculation of F(a_i) */
		fai = compute_cost(states_cur, control_cur, numsteps, T);

		/* If one of the below conditions is satisfied, bracketing phase ends */
		if (fai > f0 + a_cur*rho*df0 || fai >= fai_1){
			a1 = a_prev;
			a2 = a_cur;
			//fprintf(stderr, "-- condition 1: End of bracketing phase with a1: %lf and a2: %lf \n", a1, a2);
			break;
		}
		else{
			/* Otherwise */
			/* Compute current costates */
			compute_costates(states_cur, control_cur, numsteps, T, costates_cur);
			/* Compute current gradient */
			compute_reduced_gradient(states_cur, costates_cur, control_cur, numsteps, T, reduced_g_cur);
			/* Compute F'(a_i) */
			dfai = inner_product(reduced_g_cur, s_dir, numsteps, NU);

			/* If the below condition is satisfied, line search phase ends */
			if( fabs(dfai) <= -sigma*df0 ){
				a1 = a_cur;
				a2 = a_cur;
				//fprintf(stderr, "-- condition 2: End of line search phase with alpha: %lf \n", a_cur);
				return a_cur;
			}
			else if (dfai >= 0.0){
				a1 = a_prev;
				a2 = a_cur;
				//fprintf(stderr, "-- condition 3: End of bracketing phase with a1: %lf and a2: %lf \n", a1, a2);
				break;
			}
			a_temp = a_cur;
			a_cur = a_cur + tau * (a_cur - a_prev);
			a_prev = a_temp;

			//a1 = a_prev;
			//a2 = a_cur;
		}
	}

	/* Sectioning */
	/*            */
	int count = 0;
	while(true){
		if (fabs(a1-a2) < 10e-12 || count == 10.0*lin_it){
			if (count == 10.0*lin_it){
				fprintf(stderr, "Bracketing upper bound reached \n");
			}

			return ((a1+a2)/2.0);
		}

		/* calculation for a1 */
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				control_a1[i][j] = control[i][j] + a1*s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(control_a1, c_min, c_max, numsteps, NU);
		}

		compute_states(states, control_a1, numsteps, T, states_a1);
		fa1 = compute_cost(states_a1, control_a1, numsteps, T);
		compute_costates(states_a1, control_a1, numsteps, T, costates_a1);
		compute_reduced_gradient(states_a1, costates_a1, control_a1, numsteps, T, reduced_g_a1);
		dfa1 = inner_product(reduced_g_a1, s_dir, numsteps, NU);

		/* calculation for a2 */
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				control_a2[i][j] = control[i][j] + a2*s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(control_a2, c_min, c_max, numsteps, NU);
		}

		compute_states(states, control_a2, numsteps, T, states_a2);
		fa2 = compute_cost(states_a2, control_a2, numsteps, T);
		compute_costates(states_a2, control_a2, numsteps, T, costates_a2);
		compute_reduced_gradient(states_a2, costates_a2, control_a2, numsteps, T, reduced_g_a2);
		dfa2 = inner_product(reduced_g_a2, s_dir, numsteps, NU);

		/*  Interpolation Methods for a1 and a2 */
		//alpha = cubic_interpolation(a1, a2, fa1, fa2, dfa1, dfa2);
		alpha = quadratic_interpolation(a1, a2, dfa1, dfa2);

		/* calculation for alpha */
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				control_alpha[i][j] = control[i][j] + alpha*s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(control_alpha, c_min, c_max, numsteps, NU);
		}

		compute_states(states, control_alpha, numsteps, T, states_alpha);
		fa = compute_cost(states_alpha, control_alpha, numsteps, T);

		if (fa > f0 + rho*alpha*df0 || fa >= fa1){
			a2 = alpha;
		}
		else{
			compute_costates(states_alpha, control_alpha, numsteps, T, costates_alpha);
			compute_reduced_gradient(states_alpha, costates_alpha, control_alpha, numsteps, T, reduced_g_alpha);
			dfa = inner_product(reduced_g_alpha, s_dir, numsteps, NU);

			if (fabs(dfa) <= -sigma*df0){
				return alpha;
			}
			else if (dfa > 0.0){
				a2 = alpha;
			}
			else if (dfa < 0.0){
				a1 = alpha;
			}
		}
		count += 1;
	}
}

static double **fletcher_reeves(double **reduced_g, double **reduced_g_prev, double **s_dir, int numsteps){
	int i,j;
	double beta;

	beta = inner_product(reduced_g, reduced_g, numsteps, NU) / inner_product(reduced_g_prev, reduced_g_prev, numsteps, NU);

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			s_dir[i][j] = -1.0 * reduced_g[i][j] + beta * s_dir[i][j];
		}
	}

	return s_dir;
}

static double **polak_ribiere(double **reduced_g, double **reduced_g_prev, double **s_dir, int numsteps){
	int i,j;
	double beta;

	for (i = 0; i < NU; i++){
		for (j = 0; j < numsteps; j++){
			grad_dif[i][j] = reduced_g[i][j] - reduced_g_prev[i][j];
		}
	}

	beta = inner_product(grad_dif, reduced_g, numsteps, NU) / inner_product(reduced_g_prev, reduced_g_prev, numsteps, NU);

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			s_dir[i][j] = -1.0 * reduced_g[i][j] + beta * s_dir[i][j];
		}
	}

	return s_dir;
}

static double **DFP(double **control, double **control_prev, double **reduced_g, double **reduced_g_prev, double **s_dir, int numsteps, int it) {
	int i, j, k, l;
	double **delta, **y, ***v, ***w, **zeta;

	assert((delta = (double**)calloc(sizeof(double*), NU)));
	for (k = 0; k < NU; k++) {
		assert((delta[k] = (double*)calloc(sizeof(double), numsteps)));
	}
	assert((y = (double**)calloc(sizeof(double*), NU)));
	for (k = 0; k < NU; k++) {
		assert((y[k] = (double*)calloc(sizeof(double), numsteps)));
	}
	assert((v = (double***)calloc(sizeof(double*), AP.iter)));
	for (k = 0; k < AP.iter; k++) {
		assert((v[k] = (double**)calloc(sizeof(double*), NU)));
		for (l = 0; l < NU; l++) {
			assert((v[k][l] = (double*)calloc(sizeof(double), numsteps)));
		}
	}
	assert((w = (double***)calloc(sizeof(double*), AP.iter)));
	for (k = 0; k < AP.iter; k++) {
		assert((w[k] = (double**)calloc(sizeof(double*), NU)));
		for (l = 0; l < NU; l++) {
			assert((w[k][l] = (double*)calloc(sizeof(double), numsteps)));
		}
	}
	assert((zeta = (double**)calloc(sizeof(double*), NU)));
	for (k = 0; k < NU; k++) {
		assert((zeta[k] = (double*)calloc(sizeof(double), numsteps)));
	}

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			delta[i][j] = control[i][j] - control_prev[i][j];
			y[i][j] = reduced_g[i][j] - reduced_g_prev[i][j];
		}
	}
	
	for (k = 0; k < it - 2; k++) {
		for (i = 0; i < NU; i++) {
			for (j = 0; j < numsteps; j++) {
				sum_term[i][j] += inner_product(v[k], y, numsteps, NU)*v[k][i][j] - inner_product(w[k], y, numsteps, NU)*w[k][i][j];
			}
		}
	}

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			zeta[i][j] = y[i][j] + sum_term[i][j];
		}
	}

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			v[it-1][i][j] = delta[i][j] / sqrt(inner_product(delta, y, numsteps, NU));
			w[it-1][i][j] = zeta[i][j] / sqrt(inner_product(zeta, y, numsteps, NU));
		}
	}

	for (k = 0; k < it - 1; k++) {
		for (i = 0; i < NU; i++) {
			for (j = 0; j < numsteps; j++) {
				sum_term_2[i][j] += inner_product(v[k], reduced_g, numsteps, NU)*v[k][i][j] - inner_product(w[k], reduced_g, numsteps, NU)*w[k][i][j];
			}
		}
	}

	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			s_dir[i][j] = -reduced_g[i][j] - sum_term_2[i][j];
		}
	}

	return s_dir;
}

static void optimization(const char *criterion, double **states, double **control, double **reduced_g, int numsteps, double T, int iterations, double accuracy, double a, double sigma, double rho, double tau, int lin_it, int restart_it, double *c_min, double *c_max){

	int i, j, it;
	double cost, temp_cost, alpha;
	double beta = 0.0;

	cost = compute_cost(states, control, numsteps, T);
	for(i = 0; i < NU; i++){
		for (j = 0; j < numsteps; j++)
			reduced_g_prev[i][j] = reduced_g[i][j];
	}
	
	//fprintf(stderr, "iter: %d \t cost: %lf -- gradient: %.9f \n", 0, cost, inner_product(reduced_g, reduced_g, numsteps, NU));
	for(it = 0; it < iterations; it++){
		/* Searh Direction is Steepest Descent in the 1st iteration, if the restart  */
		/* condition is satisfied or after a pre-specified no. of iterations */
		if (it == 0 || restart(s_dir, reduced_g, numsteps, NU) || it % restart_it == 0){
			/* Steepest Descent */
			for (i = 0; i < NU; i++){
				for (j = 0; j < numsteps; j++){
					s_dir[i][j] = -1.0 * reduced_g[i][j];
				}
			}
		}
		else{
			/* Calculate search direction based on the chosen algorithm */
			switch (SDIR) {
			case 1:
				/* Steepest Descent */
				for (i = 0; i < NU; i++) {
					for (j = 0; j < numsteps; j++) {
						s_dir[i][j] = -1.0 * reduced_g[i][j];
					}
				}
				break;
			case 2:
				/* Fletcher - Reeves */
				s_dir = fletcher_reeves(reduced_g, reduced_g_prev, s_dir, numsteps);
				break;
			case 3:
				/* Polak - Ribiere */
				s_dir = polak_ribiere(reduced_g, reduced_g_prev, s_dir, numsteps);
				break;
			case 4:
				/* DFP */
				DFP(control, control_pre, reduced_g, reduced_g_prev, s_dir, numsteps, it);
				break;
			default:
				printf("/* Invalid search direction method */\n");
				return;
			}
		}
		/* if the restart condition is satisfied break and start a new iteration */
		/* with Descent Steepest algorithm */
		if (it != 0 && restart(s_dir, reduced_g, numsteps, NU)){
			continue;
		}

		/* Calculate the optimum scalar step (alpha) through Line Optimization */
		alpha = line_search(a, states, control, s_dir, sigma, tau, rho, lin_it, numsteps, T);

		/* Calculate the new control */
		for (i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++){
				temp_control[i][j] = control[i][j] + alpha * s_dir[i][j];
			}
		}
		if (CONSTRAINTS){
			control_bound(temp_control, c_min, c_max, numsteps, NU);
		}

		/* Calculate the new states and cost */
		compute_states(states, temp_control, numsteps, T, temp_states);
		temp_cost = compute_cost(temp_states, temp_control, numsteps, T);
		compute_costates(temp_states, temp_control, numsteps, T, temp_costates);
		compute_reduced_gradient(temp_states, temp_costates, temp_control, numsteps, T, temp_reduced_g);

		/* Check if the convergence test is satisfied */
		if (termination(criterion, accuracy, temp_cost, cost, temp_reduced_g, numsteps) || it == iterations - 1){
			//fprintf(stderr, "termination \n");
			for (i = 0; i < NX; i++){
				for (j = 0; j < numsteps + 1; j++)
					x_opt[i][j] = temp_states[i][j];
				//memcpy(&x_opt[i], &temp_states[i], sizeof(temp_states[0]));
			}
			for (i = 0; i < NU; i++){
				for (j = 0; j < numsteps; j++)
					u_opt[i][j] = temp_control[i][j];
				//memcpy(&u_opt[i], &temp_control[i], sizeof(temp_states[0]));
			}
			cost = temp_cost;

			/* print No. of iterations and final norm of gradient */
			fprintf(stderr, "iter: %d \t cost: %lf -- gradient: %.9f \n", it+1, cost, inner_product(temp_reduced_g, temp_reduced_g, numsteps, NU));

			break;
		}
		else{
			for(i = 0; i < NU; i++){
				for (j = 0; j < numsteps; j++) {
					reduced_g_prev[i][j] = reduced_g[i][j];
					reduced_g[i][j] = temp_reduced_g[i][j];
					control_pre[i][j] = control[i][j];    /* new for bfgs */
					control[i][j] = temp_control[i][j];
				}
				//memcpy(&reduced_g_prev[i], &reduced_g[i], sizeof(reduced_g[0]));
				//memcpy(&reduced_g[i], &temp_reduced_g[i], sizeof(reduced_g[0]));
				//memcpy(&control_pre[i], &control[i], sizeof(temp_control[0]));    /* new for bfgs */
				//memcpy(&control[i], &temp_control[i], sizeof(temp_control[0]));
			}
			for(i = 0; i < NX; i++){
				for (j = 0; j < numsteps + 1; j++)
					states[i][j] = temp_states[i][j];
				//memcpy(&states[i], &temp_states[i], sizeof(temp_states[0]));
			}
			cost = temp_cost;
		}
		//fprintf(stderr, "iter: %d \t cost: %lf -- gradient: %.9f \n", it+1, cost, inner_product(reduced_g, reduced_g, numsteps, NU));
	}
}

static void allocations(int numsteps) {
	int i;

	assert((init.u = (double**)calloc(sizeof(double*), NU)));
	assert((control = (double**)calloc(sizeof(double*), NU)));
	assert((sum_term = (double**)calloc(sizeof(double*), NU)));
	assert((sum_term_2 = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g = (double**)calloc(sizeof(double*), NU)));
	assert((control_prev = (double**)calloc(sizeof(double*), NU)));
	assert((control_cur = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g_cur = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g_a1 = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g_a2 = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g_alpha = (double**)calloc(sizeof(double*), NU)));
	assert((control_a1 = (double**)calloc(sizeof(double*), NU)));
	assert((control_a2 = (double**)calloc(sizeof(double*), NU)));
	assert((control_alpha = (double**)calloc(sizeof(double*), NU)));
	assert((grad_dif = (double**)calloc(sizeof(double*), NU)));
	assert((reduced_g_prev = (double**)calloc(sizeof(double*), NU)));
	assert((control_pre = (double**)calloc(sizeof(double*), NU)));
	assert((s_dir = (double**)calloc(sizeof(double*), NU)));
	assert((temp_control = (double**)calloc(sizeof(double*), NU)));
	assert((temp_reduced_g = (double**)calloc(sizeof(double*), NU)));
	assert((u_opt = (double**)calloc(sizeof(double*), NU)));
	for (i = 0; i < NU; i++) {
		assert((init.u[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((sum_term[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((sum_term_2[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((control_prev[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control_cur[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g_cur[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g_a1[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g_a2[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g_alpha[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control_a1[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control_a2[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control_alpha[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((grad_dif[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((reduced_g_prev[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((control_pre[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((s_dir[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((temp_control[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((temp_reduced_g[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((u_opt[i] = (double*)calloc(sizeof(double), numsteps)));
	}

	assert((init.x = (double**)calloc(sizeof(double*), NX)));
	assert((states = (double**)calloc(sizeof(double*), NX)));
	assert((costates = (double**)calloc(sizeof(double*), NX)));
	assert((states_prev = (double**)calloc(sizeof(double*), NX)));
	assert((states_cur = (double**)calloc(sizeof(double*), NX)));
	assert((costates_cur = (double**)calloc(sizeof(double*), NX)));
	assert((states_a1 = (double**)calloc(sizeof(double*), NX)));
	assert((states_a2 = (double**)calloc(sizeof(double*), NX)));
	assert((states_alpha = (double**)calloc(sizeof(double*), NX)));
	assert((costates_a1 = (double**)calloc(sizeof(double*), NX)));
	assert((costates_a2 = (double**)calloc(sizeof(double*), NX)));
	assert((costates_alpha = (double**)calloc(sizeof(double*), NX)));
	assert((x_opt = (double**)calloc(sizeof(double*), NX)));
	assert((temp_states = (double**)calloc(sizeof(double*), NX)));
	assert((temp_costates = (double**)calloc(sizeof(double*), NX)));
	for (i = 0; i < NX; i++) {
		assert((init.x[i] = (double*)calloc(sizeof(double), numsteps)));
		assert((states[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((costates[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((states_prev[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((states_cur[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((costates_cur[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((states_a1[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((states_a2[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((states_alpha[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((costates_a1[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((costates_a2[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((costates_alpha[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((x_opt[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((temp_states[i] = (double*)calloc(sizeof(double), numsteps + 1)));
		assert((temp_costates[i] = (double*)calloc(sizeof(double), numsteps + 1)));
	}
}

static void free() {

	int i, j;

	for (j = 0; j < NU; j++) {
		free(init.u[j]);
		free(control[j]);
		free(sum_term[j]);
		free(sum_term_2[j]);

		free(reduced_g[j]);
		free(control_prev[j]);
		free(control_cur[j]);
		free(reduced_g_a1[j]);
		free(reduced_g_a2[j]);
		free(reduced_g_alpha[j]);
		
		free(control_a1[j]);
		free(control_a2[j]);
		free(control_alpha[j]);
		free(grad_dif[j]);
		free(control_pre[j]);
		free(s_dir[j]);
		
		free(temp_control[j]);
		free(temp_reduced_g[j]);
		free(u_opt[j]);
	}
	free(init.u);
	free(control);
	free(sum_term);
	free(sum_term_2);

	free(reduced_g);
	free(control_prev);
	free(control_cur);
	free(reduced_g_a1);
	free(reduced_g_a2);
	free(reduced_g_alpha);
	free(control_a1);
	free(control_a2);
	free(control_alpha);
	free(grad_dif);
	free(control_pre);
	free(s_dir);

	free(temp_control);
	free(temp_reduced_g);
	free(u_opt);


	for (j = 0; j < NX; j++) {
		free(init.x[j]);
		free(states[j]);
		free(costates[j]);
		free(states_prev[j]);

		free(states_cur[j]);
		free(costates_cur[j]);
		free(states_a1[j]);
		free(states_a2[j]);
		free(states_alpha[j]);
		free(costates_a1[j]);

		free(costates_a2[j]);
		free(costates_alpha[j]);
		free(x_opt[j]);
		free(temp_states[j]);
		free(temp_costates[j]);
	}
	free(init.x);
	free(states);
	free(costates);
	free(states_prev);

	free(states_cur);
	free(costates_cur);
	free(states_a1);
	free(states_a2);
	free(states_alpha);
	free(costates_a1);

	free(costates_a2);
	free(costates_alpha);
	free(x_opt);
	free(temp_states);
	free(temp_costates);	
}

static void printsol(double **states, double **control, int numsteps, double T) {
	int i, k;
	double lane;

	puts("var Data = {");
	for (i = 0; i < P.numsteps + 1; i++) {
		double thisy;

		/* FDA optimal trajectories */
		printf("\"x(%d)\":%.5f,\n", i, states[0][i]);
		thisy = ((((3.0 * P.numlanes) - 1.5) - (states[1][i])) / 3.0);
		//thisy = states[1][i];
		printf("\"y(%d)\":%.5f,\n", i, (thisy));
		printf("\"vx(%d)\":%.5f,\n", i, states[2][i]);
		printf("\"vy(%d)\":%.5f,\n", i, states[3][i]);
		printf("\"ax(%d)\":%.5f,\n", i, states[4][i]);
		printf("\"ux(%d)\":%.5f,\n", i, control[0][i]);
		printf("\"uy(%d)\":%.5f,\n", i, control[1][i]);
		//printf("\"init_uy(%d)\":%.5f,\n", i, init.u[1][i]);
		//printf("\"init_ax(%d)\":%.5f,\n", i, P.ux0[i]);
		//printf("\"init_ux(%d)\":%.5f,\n", i, init.u[0][i]);
	}

	for (i = 0; i < P.n; i++) {
		printf("\"obst_id(%d)\":%d,\n", i, P.obst_id[i]);
		printf("\"obst_connected(%d)\":%d,\n", i, P.obst_connected[i][0]);
		printf("\"obst_vd(%d)\":%.5f,\n", i, P.obst_vd[i]);
		for (k = 0; k < P.numsteps + 1; k++) {		
			printf("\"obst_x(%d,%d)\":%.5f,\n", i, k, obs.x[i][0][k]);
			lane = (((3.0 * P.numlanes) - 1.5) - obs.x[i][1][k]) / 3.0;
			//lane = obs.x[i][1][k];
			printf("\"obst_y(%d,%d)\":%.5f,\n", i, k, lane);
			if (P.obst_connected[i][0] == 1.0)
				printf("\"obst_v(%d,%d)\":%.5f,\n", i, k, obs.v[i][k]);
			else
				printf("\"obst_v(%d,%d)\":%.5f,\n", i, k, obs.v[i][0]);
			printf("\"obst_ellipse_x(%d,%d)\":%.5f,\n", i, k, r[i][k]);
			printf("\"obst_ellipse_y(%d,%d)\":%.5f,\n", i, k, r2);
			printf("\"obst_ellipse_shift(%d,%d)\":%.5f,\n", i, k, shift[i][k]);
			puts("");
		}
		puts("");
	}

	printf("\"n\":%d,\n", P.n);
	printf("\"k\":%d,\n", P.numsteps);
	printf("\"numlanes\":%d,\n", P.numlanes);
	printf("\"vd\":%f,\n", P.vd);
	printf("\"Step\":%f,\n", P.T);
	printf("\"C\":%f,\n", P.C);
	printf("\"safety\":%f,\n", P.safety);
	printf("\"id\":%d,\n", P.id);
	printf("};\n\n");
}

static void read_P(void) {

	int i;
	int bufSize = 512;
	char buf[512];
	const char *line;
	double dval;
	int ival, ival2, ival3;

	assert((P.obst_x = (double**)calloc(sizeof(double), bufSize)));
	assert((P.obst_v = (double**)calloc(sizeof(double), bufSize)));
	assert((P.obst_y = (double**)calloc(sizeof(double), bufSize)));
	assert((P.obst_connected = (int**)calloc(sizeof(int*), bufSize)));
	assert((P.obst_ax = (double**)calloc(sizeof(double), bufSize)));
	for (i = 0; i < bufSize; i++) {
		assert((P.obst_x[i] = (double*)calloc(sizeof(double), bufSize)));
		assert((P.obst_v[i] = (double*)calloc(sizeof(double), bufSize)));
		assert((P.obst_y[i] = (double*)calloc(sizeof(double), bufSize)));
		assert((P.obst_connected[i] = (int*)calloc(sizeof(int), bufSize)));
		assert((P.obst_ax[i] = (double*)calloc(sizeof(double), bufSize)));
	}
	assert((P.obst_planTime = (double*)calloc(sizeof(double), bufSize)));
	assert((P.obst_vd = (double*)calloc(sizeof(double), bufSize)));
	assert((P.obst_id = (int*)calloc(sizeof(int), bufSize)));		
	assert((P.ux0 = (double*)calloc(sizeof(double), bufSize)));
	assert((P.uy0 = (double*)calloc(sizeof(double), bufSize)));
	assert((P.xk = (double*)calloc(sizeof(double), bufSize)));
	assert((P.yk = (double*)calloc(sizeof(double), bufSize)));
	assert((P.vk = (double*)calloc(sizeof(double), bufSize)));
	assert((P.ax_init = (double*)calloc(sizeof(double), bufSize)));
	
	
	
		while (fgets(buf, sizeof(buf), stdin))
		{
			//fprintf(stderr, "parsing: '%s'\n", buf);
			if (buf[0] == '\n')
				break;
			//fprintf(stderr, "parsing: '%s'\n", line);
			if (sscanf_s(buf, "\"vd\":%lf", &dval) == 1) {
				P.vd = dval;
			}
			if (sscanf_s(buf, "\"numlanes\":%d", &ival) == 1)
				P.numlanes = ival;
			if (sscanf_s(buf, "\"numsteps\":%d", &ival) == 1)
				P.numsteps = ival;
			if (sscanf_s(buf, "\"C\":%lf", &dval) == 1)
				P.C = dval;
			if (sscanf_s(buf, "\"T\":%lf", &dval) == 1)
				P.T = dval;
			if (sscanf_s(buf, "\"safety\":%lf", &dval) == 1) {
				P.safety = dval * smu;
			}
			if (sscanf_s(buf, "\"planTime\":%lf", &dval) == 1) {
				P.planTime = dval;
			}
			if (sscanf_s(buf, "\"id\":%d", &ival) == 1) {
				P.id = ival;
			}
			if (sscanf_s(buf, "\"crash\":%d", &ival) == 1) {
				P.init_crash = ival;
			}

			/* read initial states */
			if (sscanf_s(buf, "\"x(%d)\":%lf", &ival, &dval) == 2) {
				if (ival == 0)
					P.x0 = dval;
				P.xk[ival] = dval;
			}
			if (sscanf_s(buf, "\"y(%d)\":%lf", &ival, &dval) == 2) {
				if (ival == 0)
					P.y0 = dval;
				P.yk[ival] = dval;
			}
			if (sscanf_s(buf, "\"v(0)\":%lf", &dval) == 1) {
				if (ival == 0)
					P.v0 = dval;
				P.vk[ival] = dval;
			}
			if (sscanf_s(buf, "\"a(0)\":%lf", &dval) == 1)
				P.a0 = dval;
			/* read control (jerk) */
			if (sscanf_s(buf, "\"ux(%d)\":%lf", &ival, &dval) == 2) {
				//P.ux0[ival] = (ival == 0) ? dval : dval - P.ux0[ival - 1];
				P.ux0[ival] = dval;
			}
			if (sscanf_s(buf, "\"uy(%d)\":%lf", &ival, &dval) == 2) {
				//P.uy0[ival] = (ival == 0) ? 0 : dval - P.uy0[ival - 1];
				P.uy0[ival] = dval;
			}

			/* read obstacles */
			if (sscanf_s(buf, "\"obst_x(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_x[ival][ival2] = dval;
				P.obst_x_len = MAX(P.obst_x_len, ival + 1);
			}
			if (sscanf_s(buf, "\"obst_y(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_y[ival][ival2] = dval;
				P.obst_y_len = MAX(P.obst_y_len, ival + 1);
			}
			if (sscanf_s(buf, "\"obst_v(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_v[ival][ival2] = dval;
				P.obst_v_len = MAX(P.obst_v_len, ival + 1);
			}
			if (sscanf_s(buf, "\"obst_connected(%d,%d)\":%d", &ival, &ival2, &ival3) == 3) {
				P.obst_connected[ival][ival2] = ival3;
				//P.obst_v_len = MAX(P.obst_v_len, ival + 1);
			}
			if (sscanf_s(buf, "\"obst_planTime(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_planTime[ival] = dval;
			}
			if (sscanf_s(buf, "\"obst_vd(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_vd[ival] = dval;
			}
			if (sscanf_s(buf, "\"obst_id(%d,%d)\":%d", &ival, &ival2, &ival3) == 3) {
				P.obst_id[ival] = ival3;
			}
			if (sscanf_s(buf, "\"obst_ax(%d,%d)\":%lf", &ival, &ival2, &dval) == 3) {
				P.obst_ax[ival][ival2] = dval;
				//P.obst_v_len = MAX(P.obst_v_len, ival + 1);
			}
			if (sscanf_s(buf, "\"a(0)\":%lf", &dval) == 1)
				P.a0 = dval;

			memset(buf, 0, sizeof(buf));
		}
	
	if (!(P.obst_x_len == P.obst_y_len && P.obst_x_len == P.obst_v_len)) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	P.n = P.obst_x_len;
}

int main(){
	int i;
	
	clock_t start, start_FDA, end, end_FDA;
	double cpu_time_used, cpu_time_used_FDA;
	start = clock();
	read_P();

	allocations(P.numsteps);
	obst_trajectories(P.numsteps, P.T);

	start_FDA = clock();

	init.x[0][0] = P.x0;
	init.x[1][0] = ((3.0 * (P.y0 - 1)) + 1.5);
	init.x[2][0] = P.v0;
	init.x[3][0] = 0.0;
	init.x[4][0] = P.a0;

	P.init_crash = 1;
	for (i = 0; i < P.numsteps; i++) {
		init.u[0][i] = 0.0;
		init.u[1][i] = 0.0;
		if (DP) {
			if (P.init_crash == 1)
				init.u[0][i] = (i == 0) ? -6.0 : 0.0;	//CHECK
			if (i % 4 == 0)
				init.u[1][i] = -P.uy0[i] / P.T;
				//init.u[0][i] = (i == 0) ? 0.0 : -(P.uy0[i] - P.uy0[i - 1]) / P.T;	//CHECK
			
		}
	}

	compute_states(init.x, init.u, P.numsteps, P.T, states);
	compute_costates(states, init.u, P.numsteps, P.T, costates);
	compute_reduced_gradient(states, costates, init.u, P.numsteps, P.T, reduced_g);

	/* Choose Optimization Method */
	optimization(criterion, states, init.u, reduced_g, P.numsteps, P.T, AP.iter, AP.accur, AP.a, AP.sigma, AP.rho, AP.tau, AP.lin_it, AP.restart_it, c_min, c_max);
	
	end = clock();
	end_FDA = clock();
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	cpu_time_used_FDA = ((double)(end_FDA - start_FDA)) / CLOCKS_PER_SEC;

	/* Print Solution */
	printsol(x_opt, u_opt, P.numsteps, P.T);
	fprintf(stderr, "\n");
	fprintf(stderr, "CPU Time (Full Problem): %.4f sec. \n", cpu_time_used);
	fprintf(stderr, "CPU Time (FDA): %.4f sec. \n", cpu_time_used_FDA);

	free();
	fprintf(stderr, "OK\n");

	return 0; 
}
