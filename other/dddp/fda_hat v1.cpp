#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>

//#define sscanf sscanf_s

#define MAX(a,b) (a > b)?(a):(b)
#define MIN(a,b) (a < b)?(a):(b)
// #define MYSIGN(a) (a > 0.0)?(1.0):(-1.0)
#define NU 2 // no. of controls per vehicle
#define NX 4 // no. of states per vehicle

#define Power(a, b) pow(a, b)
#define Tanh(a) tanh(a)
#define Sech(a) (1.0/cosh(a))
#define Sqrt(a) sqrt(a)

#define SDIR 5

#define P1MULT 3
#define AMULT 0.5
#define BMULT 0.5
#define EPSILONW (1e-1)

double LATK1 = 0.125;
double LATK2 = 2.0 * sqrt(LATK1) - LATK1 * (0.25 / 2.0);

double LONK3 = 4.0;
double BETA = 0.03;

//                        {  UX,  UY, VDX, VDY, ELL, VCP, ATTRAC, PULL };
static double penalty[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

double c_min[] = { 0.0, 0.0 };
double c_max[] = { 1.0, 1.0 };

/* ellipse parameters */
int p1 = 2, p2 = 2, p3 = 2;

static double *grad_norm;
static double ***states_all_veh, ***control_all_veh;
/* debugging */
int ls_no_iter = 0, brac_it = 0, sect_it = 0;

struct vehicle {
	double **states, **costates;
	double **control;
	double **reduced_g;
	double **reduced_g_prev;
	double **temp_states, **temp_costates, **temp_control, **temp_reduced_g;
	
	double **xLim, **vxLim, **yLim, **vyLim;
	double *uxLB, *uxUB, *uyLB, *uyUB;

	/* rprop */
	double **delta_u;

	/* */
	int isLeader;
};

struct vehicle *veh;

static const char *criterion = "gradient";
static struct {
	int iter;
	double accur;
	double a;
	double sigma, rho, sigma2, tau;
	int lin_it, restart_it;

	/* rprop */
	double eta_minus, eta_plus;
	double *delta_u_0;
	double delta_u_min, delta_u_max;
}AP = { 0 };

static struct {
	double *vdx, *vdy, T, *ego_l, *ego_w;
	double *x0, *y0, *vx0, *vy0, *ax0, *ay0;
	int *isLeader;
	int leaderID;

	int x_len, y_len, vx_len, vy_len, ego_l_len, ego_w_len, ax_len, ay_len;

	int numsteps, n;
	double safety1, safety2, roadwidth, latmar;

	double cmux, cmuy;
	double ux_min, ux_max, uy_min, uy_max;

	double *obst_x, *obst_y, *obst_vx, *obst_vy, *obst_l, *obst_w;
	int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len, obst_l_len, obst_w_len, obst_ax_len, obst_ay_len;
	int n_obs;

}P = { 0 };

static struct {
	bool apply;
}r = { 0 };

static struct {
	double ***x;
	double ***v;
}obs = { 0 };

static struct {

	double len, width;
	int p1, p2;

}at = { 0 };


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

static void inner_product_1(double **A, double **B, int numsteps, int nu, double *result) {
	
	for (int i = 0; i < nu; i++){
		for (int j = 0; j < numsteps; j++){
			result[i] += A[i][j] * B[i][j];
		}
	}
}

static double sign(double a) {
	if (a > 0.0){
		return 1.0;
	}
	else if (a < 0.0){
		return -1.0;
	}else{
		return 0.0;
	}
}

static double quadratic_interpolation(double a, double b, double dfa, double dfb) {
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
		//fprintf(stderr, "negative sqrt in cubic interpolation, apply bisection \n");
		alpha = (a + b) / 2.0;
	}
	return alpha;
}

static bool termination(struct vehicle *veh, const char *termination, double accuracy, double cost_cur, double cost_prev, int numsteps) {
	double temp = 0.0;

	if (strncmp (termination, "relative", strlen(termination)) == 0){
		temp = (cost_cur - cost_prev) / (cost_cur + 0.001);
		fprintf(stderr, "%f, %f, %f, %f\n", cost_cur, cost_prev, temp, accuracy);
		if (fabs(temp) < accuracy){
			fprintf(stderr, "returned true\n" );
			return true;
		}
	}
	else if (strncmp (termination, "gradient", strlen(termination)) == 0){
		for (int idVeh = 0; idVeh < P.n; idVeh++)
			temp += inner_product(veh[idVeh].temp_reduced_g, veh[idVeh].temp_reduced_g, numsteps, NU);
		if (temp < accuracy){
			return true;
		}
	}
	return false;
}

static bool termination_1(struct vehicle *veh, const char *termination, double accuracy, double cost_cur, double cost_prev, int numsteps, double *grad_norm, int it) {
	double temp;

	if (strncmp (termination, "relative", strlen(termination)) == 0){
		temp = (cost_cur - cost_prev) / (cost_cur + 0.001);
		if (fabs(cost_cur - cost_prev) < 1e-12){
			return true;
		}
	}
	else if (strncmp (termination, "gradient", strlen(termination)) == 0) {
		double grad_norm_prev[] = { grad_norm[0], grad_norm[1] };
		
		grad_norm[0] = grad_norm[1] = 0.0;
		for (int idVeh = 0; idVeh < P.n; idVeh++) {
			inner_product_1(veh[idVeh].temp_reduced_g, veh[idVeh].temp_reduced_g, numsteps, NU, grad_norm);
		}
		
	// 	if ((grad_norm[0] < accuracy && grad_norm[1] < accuracy) ||
	// 	(it > 2 && ((fabs(grad_norm[0]-grad_norm_prev[0]) < 1e-4 && fabs(grad_norm[1]-grad_norm_prev[1]) < 1e-4) || fabs(cost_cur - cost_prev) < 1e-9))
	// 	|| (it > 25 && fabs(cost_cur - cost_prev) < 1e-4)) {
	// 		return true;
	// 	}
	// }
	// if ((grad_norm[0] / P.n < accuracy && grad_norm[1] / P.n < accuracy)) {
	if ((grad_norm[0] < accuracy && grad_norm[1] < accuracy)) {
			return true;
		}
	}

	return false;
}

static void obst_trajectories(int numsteps, double T) {
	for (int i = 0; i < P.n_obs; i++) {
		obs.x[i][0][0] = P.obst_x[i];
		obs.x[i][1][0] = P.obst_y[i];
		obs.v[i][0][0] = P.obst_vx[i];
		obs.v[i][1][0] = P.obst_vy[i];

	}
	for (int i = 0; i < P.n_obs; i++) {
		// if (P.obst_connected[i] == 1) {
		// 	for (int j = 0; j < numsteps; j++) {
		// 		obs.x[i][0][j + 1] = obs.x[i][0][j] + T * obs.v[i][0][j] + 0.5*pow(T, 2) * P.obst_ax[i][j];
		// 		obs.x[i][1][j + 1] = obs.x[i][1][j] + T * obs.v[i][1][j] + 0.5*pow(T, 2) * P.obst_ay[i][j];
		// 		obs.v[i][0][j + 1] = obs.v[i][0][j] + T * P.obst_ax[i][j];
		// 		obs.v[i][1][j + 1] = obs.v[i][1][j] + T * P.obst_ay[i][j];
		// 	}
		// }
		// else {
		for (int j = 0; j < numsteps; j++) {
			obs.x[i][0][j + 1] = obs.x[i][0][j] + T * obs.v[i][0][0];
			//obs.x[i][1][j + 1] = obs.x[i][1][j] + T * obs.v[i][1][0];
			obs.x[i][1][j + 1] = obs.x[i][1][j];
		}
		// }
	}
}

static void compute_states(struct vehicle* veh, double **x_init, double **control, int numsteps, double T, double **states_calc, int idVeh) {
	int i;

	for(i = 0; i < NX; i++) {
		states_calc[i][0] = x_init[i][0];
	}

	double uxTrans, uyTrans, yLB, yUB, vyLB, vyUB;
	for (i = 0; i < numsteps; i++) {
		
		yLB = veh[idVeh].yLim[0][i]; yUB = veh[idVeh].yLim[1][i]; vyLB = veh[idVeh].vyLim[0][i]; vyUB = veh[idVeh].vyLim[1][i];

		veh[idVeh].uxLB[i] = -LONK3 * states_calc[2][i];
		veh[idVeh].uxUB[i] = P.ux_max;

		veh[idVeh].uyLB[i] = -LATK1 * (states_calc[1][i] - yLB) - LATK2 * (states_calc[3][i] - vyLB);
		veh[idVeh].uyUB[i] = -LATK1 * (states_calc[1][i] - yUB) - LATK2 * (states_calc[3][i] - vyUB);

		uxTrans = (1.0 - control[0][i]) * veh[idVeh].uxLB[i] + control[0][i] * veh[idVeh].uxUB[i];
		uyTrans = (1.0 - control[1][i]) * veh[idVeh].uyLB[i] + control[1][i] * veh[idVeh].uyUB[i];
		
		states_calc[0][i + 1] = states_calc[0][i] + states_calc[2][i]*T + 0.5*uxTrans*pow(T, 2);
		states_calc[1][i + 1] = states_calc[1][i] + states_calc[3][i]*T + 0.5*uyTrans*pow(T, 2);
		states_calc[2][i + 1] = states_calc[2][i] + uxTrans*T;
		states_calc[3][i + 1] = states_calc[3][i] + uyTrans*T;
	}
}

static void compute_costates(struct vehicle* veh, double ***states, double ***control, int numsteps, double T, double **costates_calc, int idVeh) {
	double dvdx, dvdy, dvcpdvx, dvcpdvy;
	double x, y, vx, vy, ux, uy;
	
	double duxdx, duxdvx, duydy, duydvy;

	double l_m = P.ego_w[idVeh] / 2.0;
	double y0 = P.roadwidth - l_m;

	for (int i = 0; i < NX; i++) {
		costates_calc[i][numsteps] = 0.0;	/* all costates equal to zero because we dont have terminal conditions */
	}
	
	int i = numsteps - 1;
	while (i >= 0) {
		x = states[idVeh][0][i];
		y = states[idVeh][1][i];
		vx = states[idVeh][2][i];
		vy = states[idVeh][3][i];
		
		ux = control[idVeh][0][i];
		uy = control[idVeh][1][i];

		double at_costate_term_1 = 0.0, at_costate_term_2 = 0.0;
		double costate_term_1 = 0.0, costate_term_2 = 0.0, costate_term_3 = 0.0, costate_term_4 = 0.0;
		/* collision avoidance term derivatice w.r.t x, y, v */
		for (int j = idVeh + 1; j < P.n; j++) {
			double obsvx = states[j][2][i], obsvy = states[j][3][i], obsx = states[j][0][i], obsy = states[j][1][i];
			// double x = states[idVeh][0][i], y = states[idVeh][1][i], vx = states[idVeh][2][i], vy = states[idVeh][3][i];
			double safety1 = P.safety1, safety2 = P.safety2; 
			double L = 0.5*(P.ego_l[idVeh] + P.ego_l[j]), W = 0.5*(P.ego_w[idVeh] + P.ego_w[j]);
			
			costate_term_1 += -((Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) / (L + safety1 * (obsvx + vx))) -
				(1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) / (AMULT * (L + safety1 * (obsvx + vx)));

			costate_term_2 += -(Power(2, p2) * p2 * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2) *
				(-((safety2 * (-obsy + y) * (-((-obsvy + vy) * Power(Sech(obsy - y), 2)) - (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2)) + 1 / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))))) -
				1. * Power(2, p2) * p2 * p3 * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3) *
				(-((safety2 * (-obsy + y) * (-((-obsvy + vy) * Power(Sech(obsy - y), 2)) - (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2))) +
				1 / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))));

			costate_term_3 += -(Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) * (safety1 / (2. * (L + safety1 * (obsvx + vx))) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / Power(L + safety1 * (obsvx + vx), 2)) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) -
				1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				(safety1 / (2. * AMULT * (L + safety1 * (obsvx + vx))) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / (AMULT * Power(L + safety1 * (obsvx + vx), 2))) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3);

			costate_term_4 += (Power(2, p2) * p2 * safety2 * (-obsy + y) * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				(Tanh(obsy - y) + ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2)) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2) +
				(1. * Power(2, p2) * p2 * p3 * safety2 * (-obsy + y) * (Tanh(obsy - y) + ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2));
		}

		for (int j = idVeh - 1; j >= 0; j--) {
			double vx = states[j][2][i], vy = states[j][3][i], x = states[j][0][i], y = states[j][1][i];
			double obsx = states[idVeh][0][i], obsy = states[idVeh][1][i], obsvx = states[idVeh][2][i], obsvy = states[idVeh][3][i];
			double safety1 = P.safety1, safety2 = P.safety2; 
			double L = 0.5*(P.ego_l[idVeh] + P.ego_l[j]), W = 0.5*(P.ego_w[idVeh] + P.ego_w[j]);

			costate_term_1 += (Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) / (L + safety1 * (obsvx + vx)) +
				(1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) / (AMULT * (L + safety1 * (obsvx + vx)));

			costate_term_2 += -(Power(2, p2) * p2 * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2) *
				(-((safety2 * (-obsy + y) * ((-obsvy + vy) * Power(Sech(obsy - y), 2) + (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2)) - 1 / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))))) -
				1. * Power(2, p2) * p2 * p3 * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3) *
				(-((safety2 * (-obsy + y) * ((-obsvy + vy) * Power(Sech(obsy - y), 2) + (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2))) -
				1 / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))));

			costate_term_3 += -(Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) * (-0.5 * safety1 / (L + safety1 * (obsvx + vx)) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / Power(L + safety1 * (obsvx + vx), 2)) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) -
				1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				(-0.5 * safety1 / (AMULT * (L + safety1 * (obsvx + vx))) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / (AMULT * Power(L + safety1 * (obsvx + vx), 2))) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3);

			costate_term_4 += (Power(2, p2) * p2 * safety2 * (-obsy + y) * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				(-Tanh(obsy - y) - ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2)) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2) +
				(1. * Power(2, p2) * p2 * p3 * safety2 * (-obsy + y) * (-Tanh(obsy - y) - ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2));
		}

		double refx = states[P.leaderID][0][i], refy = states[P.leaderID][1][i];
		double len = at.len, width = at.width, vehLength = P.ego_w[P.leaderID];
		if (P.leaderID != idVeh) {
			at_costate_term_1 = (at.p1*Power((len - refx - vehLength + x)/len, -1 + at.p1)) / (len*Power(1 + Power((len - refx - vehLength + x)/len, at.p1) + Power((-refy + y)/width, at.p2), 2));
			at_costate_term_2 = (at.p2*Power((-refy + y)/width, -1 + at.p2)) / (width*Power(1 + Power((len - refx - vehLength + x)/len, at.p1) + Power((-refy + y)/width, at.p2), 2));
		}

		for (int ob = 0; ob < P.n_obs; ob++) {
			double obsvx = obs.v[ob][0][i], obsvy = obs.v[ob][1][i], obsx = obs.x[ob][0][i], obsy = obs.x[ob][1][i];
			double x = states[idVeh][0][i], y = states[idVeh][1][i], vx = states[idVeh][2][i], vy = states[idVeh][3][i];
			double L = 0.5*(P.ego_l[idVeh] + P.obst_l[ob]), W = 0.5*(P.ego_w[idVeh] + P.obst_w[ob]);
			double safety1 = P.safety1, safety2 = P.safety2; 

			costate_term_1 += -((Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) / (L + safety1 * (obsvx + vx))) -
				(1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) / (AMULT * (L + safety1 * (obsvx + vx)));

			costate_term_2 += -(Power(2, p2) * p2 * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2) *
				(-((safety2 * (-obsy + y) * (-((-obsvy + vy) * Power(Sech(obsy - y), 2)) - (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2)) + 1 / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))))) -
				1. * Power(2, p2) * p2 * p3 * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3) *
				(-((safety2 * (-obsy + y) * (-((-obsvy + vy) * Power(Sech(obsy - y), 2)) - (Power(-obsvy + vy, 2) * Power(Sech(obsy - y), 2) * Tanh(obsy - y)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2))) +
				1 / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))));

			costate_term_3 += -(Power(2, p1 * P1MULT) * p1 * P1MULT * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), -1 + p1 * P1MULT) * (safety1 / (2. * (L + safety1 * (obsvx + vx))) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / Power(L + safety1 * (obsvx + vx), 2)) *
				Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2)) -
				1. * Power(2, p1) * p1 * p3 * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), -1 + p1) *
				(safety1 / (2. * AMULT * (L + safety1 * (obsvx + vx))) - (safety1 * (-obsx + (safety1 * (-obsvx + vx)) / 2. + x)) / (AMULT * Power(L + safety1 * (obsvx + vx), 2))) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3);

			costate_term_4 += (Power(2, p2) * p2 * safety2 * (-obsy + y) * Power(Sech(Power(2, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
				Power(2, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2)), 2) *
				(Tanh(obsy - y) + ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), -1 + p2)) /
				Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2) +
				(1. * Power(2, p2) * p2 * p3 * safety2 * (-obsy + y) * (Tanh(obsy - y) + ((-obsvy + vy) * Power(Tanh(obsy - y), 2)) / Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))) *
				Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), -1 + p2) *
				Power(1 + Power(2, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2. + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), -1 - p3)) /
				(BMULT * Power(W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))), 2));
		}

		dvdx = penalty[2] * (states[idVeh][2][i] - P.vdx[idVeh]) / pow(P.vdx[idVeh], 2);
		dvdy = penalty[3] * (states[idVeh][3][i] - P.vdy[idVeh]);
		
		dvcpdvx = dvcpdvy = 0.0;
		if (states[idVeh][3][i] > 0 && BETA*states[idVeh][2][i] < states[idVeh][3][i]) {
			dvcpdvx = penalty[5] * (BETA*states[idVeh][2][i] - states[idVeh][3][i]) * BETA;
			dvcpdvy = -penalty[5] * (BETA*states[idVeh][2][i] - states[idVeh][3][i]);
		}
		if (states[idVeh][3][i] < 0 && BETA*states[idVeh][2][i] < -states[idVeh][3][i]) {
			dvcpdvx = penalty[5] * (BETA*states[idVeh][2][i] + states[idVeh][3][i]) * BETA;
			dvcpdvy = penalty[5] * (BETA*states[idVeh][2][i] + states[idVeh][3][i]);
		}

		double dudx, dudvx, dudy, dudvy;
		double yLB = veh[idVeh].yLim[0][i], yUB = veh[idVeh].yLim[1][i], vyLB = veh[idVeh].vyLim[0][i], vyUB = veh[idVeh].vyLim[1][i];

		dudx = 0.0;
		dudy = 1. * penalty[1] * (-(LATK1 * (1 - uy)) - LATK1 * uy) * ((1 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) +
			uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB)));
		dudvx = -1. * LONK3 * penalty[0] * (1 - ux) * (ux * P.ux_max - LONK3 * (1 - ux) * vx);
		dudvy = 1. * penalty[1] * (-(LATK2 * (1 - uy)) - LATK2 * uy) * ((1 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) +
			uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB)));


		double a, b;
		// double shift = len - vehLength;
		double xL = refx - (2*len), xU = refx, yL = refy - width, yU = refy + width;
		double DpullxL = 0.0, DpullxU = 0.0, DpullyL = 0.0, DpullyU = 0.0;
		
		if (P.leaderID != idVeh) {
			DpullxL = penalty[7] * (1.*(-1 - (0. - x + xL)/Sqrt(0.01 + Power(0. - x + xL,2)))*(0. - x + xL + Sqrt(0.01 + Power(0. - x + xL,2))));
			DpullxU = penalty[7] * (1.*(1 + (0. + x - xU)/Sqrt(0.01 + Power(0. + x - xU,2)))*(0. + x + Sqrt(0.01 + Power(0. + x - xU,2)) - xU));
			DpullyL = penalty[7] * (1.*(-1 - (0. - y + yL)/Sqrt(0.01 + Power(0. - y + yL,2)))*(0. - y + yL + Sqrt(0.01 + Power(0. - y + yL,2))));
			DpullyU = penalty[7] * (1.*(1 + (0. + y - yU)/Sqrt(0.01 + Power(0. + y - yU,2)))*(0. + y + Sqrt(0.01 + Power(0. + y - yU,2)) - yU));
		}

		/* costates calculation */
		costates_calc[0][i] = costates_calc[0][i + 1] + DpullxL + DpullxU
			+ penalty[4] * costate_term_1 + penalty[6] * at_costate_term_1;
		costates_calc[1][i] = penalty[4] * costate_term_2 + dudy + penalty[6] * at_costate_term_2 + DpullyL + DpullyU
			+ costates_calc[3][i + 1] * T * (-(LATK1 * (1 - uy)) - LATK1 * uy) + costates_calc[1][i + 1] * (1 + 0.5 * Power(T, 2) * (-(LATK1 * (1 - uy)) - LATK1 * uy));
		costates_calc[2][i] = dvdx + penalty[4] * costate_term_3 + dvcpdvx + dudvx
			+ costates_calc[2][i + 1] * (1 - LONK3 * T * (1 - ux)) + costates_calc[0][i + 1] * (T - 0.5 * LONK3 * Power(T, 2) * (1 - ux));
		costates_calc[3][i] = dvdy + penalty[4] * costate_term_4 + dvcpdvy + dudvy
			+ costates_calc[3][i + 1] * (1 + T * (-(LATK2 * (1 - uy)) - LATK2 * uy)) + costates_calc[1][i + 1] * (T + 0.5 * Power(T, 2) * (-(LATK2 * (1 - uy)) - LATK2 * uy));

		i = i - 1;
	}
}

static double compute_cost(struct vehicle* veh, double ***states, double ***control, int numsteps, double T, int idVeh) {
	int i, j;
	double ux, uy, vdx, vdy, ellipse_cost, vcp, attraction_cost;
	double uxTrans, uyTrans;
	double cost = 0.0;

	double l_m = P.ego_w[idVeh] / 2.0;
	double y0 = P.roadwidth - l_m;

	double x, y, vx, vy;

	for (i = 0; i < numsteps; i++){
		ellipse_cost = 0.0;
		attraction_cost = 0.0;
		x = states[idVeh][0][i]; y = states[idVeh][1][i]; vx = states[idVeh][2][i]; vy = states[idVeh][3][i];

		/* collision avoidance term */
		for (j = idVeh + 1; j < P.n; j++) {
			double obsvx = states[j][2][i], obsvy = states[j][3][i], obsx = states[j][0][i], obsy = states[j][1][i];
			double safety1 = P.safety1, safety2 = P.safety2; 
			double L = 0.5*(P.ego_l[idVeh] + P.ego_l[j]), W = 0.5*(P.ego_w[idVeh] + P.ego_w[j]);

			ellipse_cost += 1.0 + 1.0 / Power(1.0 + Power(2.0, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2.0 + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2.0, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), p3) -
				Tanh(Power(2.0, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2.0 + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
					Power(2.0, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2));
		}

		double refx = states[P.leaderID][0][i], refy = states[P.leaderID][1][i];
		double len = at.len, width = at.width, vehLength =  P.ego_w[P.leaderID];
		if (P.leaderID != idVeh)
			attraction_cost += 1 - 1/(1 + Power((len - refx - vehLength + x)/len,p1) + Power((-refy + y)/width,p2));

		for (int ob = 0; ob < P.n_obs; ob++) {
			double obsvx = obs.v[ob][0][i], obsvy = obs.v[ob][1][i], obsx = obs.x[ob][0][i], obsy = obs.x[ob][1][i];
			double x = states[idVeh][0][i], y = states[idVeh][1][i], vx = states[idVeh][2][i], vy = states[idVeh][3][i];
			double L = 0.5*(P.ego_l[idVeh] + P.obst_l[ob]), W = 0.5*(P.ego_w[idVeh] + P.obst_w[ob]);
			double safety1 = P.safety1, safety2 = P.safety2; 

			ellipse_cost += 1.0 + 1.0 / Power(1.0 + Power(2.0, p1) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2.0 + x) / (AMULT * (L + safety1 * (obsvx + vx))), p1) +
				Power(2.0, p2) * Power((-obsy + y) / (BMULT * (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2))))), p2), p3) -
				Tanh(Power(2.0, p1 * P1MULT) * Power((-obsx + (safety1 * (-obsvx + vx)) / 2.0 + x) / (L + safety1 * (obsvx + vx)), p1 * P1MULT) +
					Power(2.0, p2) * Power((-obsy + y) / (W + safety2 * ((-obsvy + vy) * Tanh(obsy - y) + Sqrt(EPSILONW + Power(-obsvy + vy, 2) * Power(Tanh(obsy - y), 2)))), p2));
		}

		ellipse_cost *= penalty[4];
		attraction_cost *= penalty[6];

		uxTrans = (1.0 - control[idVeh][0][i]) * veh[idVeh].uxLB[i] + control[idVeh][0][i] * veh[idVeh].uxUB[i];
		uyTrans = (1.0 - control[idVeh][1][i]) * veh[idVeh].uyLB[i] + control[idVeh][1][i] * veh[idVeh].uyUB[i];
		/* terms for u (longtitudinal/lateral) and desired speed */
		ux = 0.5 * penalty[0] * pow(uxTrans, 2);
		uy = 0.5 * penalty[1] * pow(uyTrans, 2);
		vdx = 0.5 * penalty[2] * pow( (states[idVeh][2][i] - P.vdx[idVeh]) / P.vdx[idVeh], 2);
		vdy = 0.5 * penalty[3] * pow(states[idVeh][3][i] - P.vdy[idVeh], 2);
		
		double a, b;
		double shift = len - vehLength;
		double xL = refx - (2*len), xU = refx, yL = refy - width, yU = refy + width;
		double pullxL = 0.0, pullxU = 0.0, pullyL = 0.0, pullyU = 0.0;

		if (P.leaderID != idVeh) {
			a = -xL + x; b = 0.0;
			double smoothPullxL = -0.5 * ( -a -b + sqrt(pow((-a+b), 2) + 0.01) );
			a = xU - x; b = 0.0;
			double smoothPullxU = -0.5 * ( -a -b + sqrt(pow((-a+b), 2) + 0.01) );
			a = -yL + y; b = 0.0;
			double smoothPullyL = -0.5 * ( -a -b + sqrt(pow((-a+b), 2) + 0.01) );
			a = yU - y; b = 0.0;
			double smoothPullyU = -0.5 * ( -a -b + sqrt(pow((-a+b), 2) + 0.01) );
			
			pullxL = 0.5 * penalty[7] * pow(smoothPullxL, 2);
			pullxU = 0.5 * penalty[7] * pow(smoothPullxU, 2);
			pullyL = 0.5 * penalty[7] * pow(smoothPullyL, 2);
			pullyU = 0.5 * penalty[7] * pow(smoothPullyU, 2);
		}
		vcp = 0.0;
		if (states[idVeh][3][i] > 0 && BETA*states[idVeh][2][i] < states[idVeh][3][i]) {
			vcp = 0.5 * penalty[5] * pow(BETA*states[idVeh][2][i] - states[idVeh][3][i], 2);
		}
		if (states[idVeh][3][i] < 0 && BETA*states[idVeh][2][i] < -states[idVeh][3][i]) {
			vcp = 0.5 * penalty[5] * pow(BETA*states[idVeh][2][i] + states[idVeh][3][i], 2);
		}

		cost += ux + uy + vdx + vdy + ellipse_cost + vcp + attraction_cost + pullxL + pullxU + pullyL + pullyU;
	}
	// fprintf(stderr, "cost: %.7f \t ux: %.7f \t uy: %.7f \t vdx: %.7f \t vdy: %.7f \t ell: %.7f \t vcp: %.7f\n", cost, ux, uy, vdx, vdy, ellipse_cost, vcp);

	return cost;
}

static void compute_reduced_gradient(struct vehicle *veh, double ***states, double ***control, int numsteps, double T, double **costates, double **reduced_g_calc, int idVeh){
	int i, j;
	double x, y, vx, vy, ux, uy;

	for (i = 0; i < numsteps; i++) {

		x = states[idVeh][0][i];
		y = states[idVeh][1][i];
		vx = states[idVeh][2][i];
		vy = states[idVeh][3][i];

		ux = control[idVeh][0][i];
		uy = control[idVeh][1][i];

		double yLB = veh[idVeh].yLim[0][i], yUB = veh[idVeh].yLim[1][i], vyLB = veh[idVeh].vyLim[0][i], vyUB = veh[idVeh].vyLim[1][i];
		// fprintf(stderr, "(%d, %d) yLB: %.4f \t yUB: %.4f \t vyLB: %.4f \t vyUB: %.4f \n", idVeh, i, yLB, yUB, vyLB, vyUB);
		reduced_g_calc[0][i] = 1.0 * penalty[0] * (P.ux_max + LONK3 * vx) * (ux * P.ux_max - LONK3 * (1.0 - ux) * vx)
			+ costates[2][i + 1] * T * (P.ux_max + LONK3 * vx) + 0.5 * costates[0][i + 1] * Power(T, 2) * (P.ux_max + LONK3 * vx);

		reduced_g_calc[1][i] = 1.0 * penalty[1] * ((1.0 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) + uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB))) *
			(LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB))
			+ costates[3][i + 1] * T * (LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB)) +
			0.5 * costates[1][i + 1] * Power(T, 2) * (LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB));
	}

	for (i = 0; i < NU; i++){
		for (j = 0; j < numsteps; j++){
			if ((control[idVeh][i][j] == c_min[i]) && (reduced_g_calc[i][j] > 0.0)) {
				reduced_g_calc[i][j] = 0.0;
			}
			else if( (control[idVeh][i][j] == c_max[i]) && (reduced_g_calc[i][j] < 0.0) ){
				reduced_g_calc[i][j] = 0.0;
			}
		}
	}
	
}

static void control_bound(double **control, double *control_min, double *control_max, int numsteps, int nu, int idVeh){
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

static void rporp(struct vehicle *veh, int numsteps, int nu) {

	for (int idVeh = 0; idVeh < P.n; idVeh++ ) {
		for(int i = 0; i < NU; i++){
			for (int j = 0; j < numsteps; j++) {
				
				if (veh[idVeh].reduced_g_prev[i][j] * veh[idVeh].reduced_g[i][j] > 0.0)
					veh[idVeh].delta_u[i][j] = MIN(AP.eta_plus * (veh[idVeh].delta_u[i][j]), AP.delta_u_max);
				else if (veh[idVeh].reduced_g_prev[i][j] * veh[idVeh].reduced_g[i][j] < 0.0)
					veh[idVeh].delta_u[i][j] = MAX(AP.eta_minus * (veh[idVeh].delta_u[i][j]), AP.delta_u_min);
				else
					veh[idVeh].delta_u[i][j] = veh[idVeh].delta_u[i][j];

				if (veh[idVeh].reduced_g[i][j] > 0.0)
					veh[idVeh].temp_control[i][j] = veh[idVeh].control[i][j] - veh[idVeh].delta_u[i][j];
				else if (veh[idVeh].reduced_g[i][j] < 0.0)
					veh[idVeh].temp_control[i][j] = veh[idVeh].control[i][j] + veh[idVeh].delta_u[i][j];
				else
					veh[idVeh].temp_control[i][j] = veh[idVeh].control[i][j];
			}
		}
	}


}

static void optimization(struct vehicle *veh, const char *criterion, int numsteps, double T, int iterations, double accuracy, double a, double sigma, double rho, double tau, int lin_it, int restart_it){

	int i, j, it, idVeh, k;
	double cost, temp_cost, alpha, cost_prev;
	double beta = 0.0;

	clock_t start_opt;
	start_opt = clock();


	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < NX; i++) {
			for (k = 0; k < P.numsteps + 1; k++) {
				states_all_veh[idVeh][i][k] = veh[idVeh].states[i][k];
			}
		}
		for (i = 0; i < NU; i++) {
			for (k = 0; k < P.numsteps; k++) {
				control_all_veh[idVeh][i][k] = veh[idVeh].control[i][k];
			}
		}
	}

	cost = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++) { cost += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh); }
	
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for(i = 0; i < NU; i++){
			for (j = 0; j < numsteps; j++) {
				veh[idVeh].reduced_g_prev[i][j] = veh[idVeh].reduced_g[i][j];
			}
		}
	}

	grad_norm[0] = grad_norm[1] = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++) { inner_product_1(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, NU, grad_norm); }
	
	fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", 0, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);
	
	FILE *fc;
	fc = fopen("cost_transform.txt", "w");
	fprintf(fc, "%d \t %.7f \t %.4f \t %.4f \t %.4f \n", 0, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC);
	
	for(it = 0; it < iterations; it++) {
		
		rporp(veh, numsteps, NU);
		
		for (idVeh = 0; idVeh < P.n; idVeh++) { control_bound(veh[idVeh].temp_control, c_min, c_max, numsteps, NU, idVeh); }

		/* Calculate the new states and cost */
		for (idVeh = 0; idVeh < P.n; idVeh++) { compute_states(veh, veh[idVeh].states, veh[idVeh].temp_control, numsteps, T, veh[idVeh].temp_states, idVeh); }

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].temp_states[i][k];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].temp_control[i][k];
				}
			}
		}

		temp_cost = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			temp_cost += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);
			compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].temp_costates, idVeh);
			compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].temp_costates, veh[idVeh].temp_reduced_g, idVeh);
		}

		/* Check if the convergence test is satisfied */
		if (termination_1(veh, criterion, accuracy, temp_cost, cost, numsteps, grad_norm, it) || it == iterations - 1){
			// fprintf(stderr, "termination \n");
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NX; i++){
					for (j = 0; j < numsteps+1; j++) {
						veh[idVeh].states[i][j] = veh[idVeh].temp_states[i][j];
					}
				}
				for (i = 0; i < NU; i++){
					for (j = 0; j < numsteps; j++) {
						veh[idVeh].control[i][j] = veh[idVeh].temp_control[i][j];
					}
				}
			}
			cost = temp_cost;

			fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", it+1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);
			fprintf(fc, "%d \t %.7f \t %.4f \t %.4f \t %.4f \n", it+1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC);
			break;
		}
		else {
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for(i = 0; i < NU; i++){
					for (j = 0; j < numsteps; j++) {
						veh[idVeh].reduced_g_prev[i][j] = veh[idVeh].reduced_g[i][j];
						veh[idVeh].reduced_g[i][j] = veh[idVeh].temp_reduced_g[i][j];
						/* added rprop */
						veh[idVeh].control[i][j] = veh[idVeh].temp_control[i][j];
					}
				}
				for(i = 0; i < NX; i++){
					for (j = 0; j < numsteps + 1; j++) {
						veh[idVeh].states[i][j] = veh[idVeh].temp_states[i][j];
					}
				}
			}
			cost_prev = cost;
			cost = temp_cost;
		}

		grad_norm[0] = grad_norm[1] = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++) { inner_product_1(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, NU, grad_norm); }

		fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", it+1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);
		fprintf(fc, "%d \t %.7f \t %.4f \t %.4f \t %.4f \n", it+1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC);
	}
}

static void allocations(int numsteps) {
	int i, j, idVeh;

	((veh = (struct vehicle*)calloc(P.n, sizeof(struct vehicle))));
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		((veh[idVeh].control = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_prev = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].temp_control = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].temp_reduced_g = (double**)calloc(NU, sizeof(double*))));

		((veh[idVeh].delta_u = (double**)calloc(NU, sizeof(double*))));
		for (i = 0; i < NU; i++) {
			((veh[idVeh].control[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_prev[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].temp_control[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].temp_reduced_g[i] = (double*)calloc(numsteps, sizeof(double))));

			((veh[idVeh].delta_u[i] = (double*)calloc(numsteps, sizeof(double))));
		}

		((veh[idVeh].uxLB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uxUB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uyLB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uyUB = (double*)calloc(numsteps, sizeof(double))));
		
		((veh[idVeh].states = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].temp_states = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].temp_costates = (double**)calloc(NX, sizeof(double*))));
		for (i = 0; i < NX; i++) {
			((veh[idVeh].states[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].temp_states[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].temp_costates[i] = (double*)calloc(numsteps + 1, sizeof(double))));
		}

		((veh[idVeh].xLim = (double**)calloc(2, sizeof(double*))));
		((veh[idVeh].vxLim = (double**)calloc(2, sizeof(double*))));
		((veh[idVeh].yLim = (double**)calloc(2, sizeof(double*))));
		((veh[idVeh].vyLim = (double**)calloc(2, sizeof(double*))));
		for (i = 0; i < 2; i++) {
			((veh[idVeh].xLim[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].vxLim[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].yLim[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].vyLim[i] = (double*)calloc(numsteps + 1, sizeof(double))));
		}
	}
	((grad_norm = (double*)calloc(NU, sizeof(double))));

	((states_all_veh = (double***)calloc(P.n, sizeof(double**))));
	for (i = 0; i < P.n; i++) {
		((states_all_veh[i] = (double**)calloc(NX, sizeof(double*))));
		for (j = 0; j < NX; j++)
			((states_all_veh[i][j] = (double*)calloc(numsteps + 1, sizeof(double))));
	}
	((control_all_veh = (double***)calloc(P.n, sizeof(double**))));
	for (i = 0; i < P.n; i++) {
		((control_all_veh[i] = (double**)calloc(NU, sizeof(double*))));
		for (j = 0; j < NU; j++)
			((control_all_veh[i][j] = (double*)calloc(numsteps, sizeof(double))));
	}

	((obs.x = (double***)calloc(sizeof(double**), P.n_obs)));
	((obs.v = (double***)calloc(sizeof(double**), P.n_obs)));
	for (int i = 0; i < P.n_obs; i++) {
		((obs.x[i] = (double**)calloc(sizeof(double*), NU)));
		((obs.v[i] = (double**)calloc(sizeof(double*), NU)));
		for (int j = 0; j < NU; j++) {
			((obs.x[i][j] = (double*)calloc(sizeof(double), numsteps + 1)));
			((obs.v[i][j] = (double*)calloc(sizeof(double), numsteps + 1)));
		}
	}
}

static void read_P(void) {

	int bufSize = 1024;
	char buf[1024];
	// const char* line;
	double dval;
	int ival, ival2;

	((P.x0 = (double*)calloc(sizeof(double), bufSize)));
	((P.y0 = (double*)calloc(sizeof(double), bufSize)));
	((P.vx0 = (double*)calloc(sizeof(double), bufSize)));
	((P.vy0 = (double*)calloc(sizeof(double), bufSize)));
	((P.ax0 = (double*)calloc(sizeof(double), bufSize)));
	((P.ay0 = (double*)calloc(sizeof(double), bufSize)));

	((P.vdx = (double*)calloc(sizeof(double), bufSize)));
	((P.vdy = (double*)calloc(sizeof(double), bufSize)));
	((P.ego_l = (double*)calloc(sizeof(double), bufSize)));
	((P.ego_w = (double*)calloc(sizeof(double), bufSize)));

	((P.obst_x = (double*)calloc(sizeof(double), bufSize)));
	((P.obst_vx = (double*)calloc(sizeof(double), bufSize)));
	((P.obst_vy = (double*)calloc(sizeof(double), bufSize)));
	((P.obst_y = (double*)calloc(sizeof(double), bufSize)));
	((P.obst_l = (double*)calloc(sizeof(double), bufSize)));
	((P.obst_w = (double*)calloc(sizeof(double), bufSize)));

	((P.isLeader = (int*)calloc(sizeof(int), bufSize)));

	((AP.delta_u_0 = (double*)calloc(sizeof(double), NU)));

	FILE* input;
	input = fopen("input-8.txt", "r");

	while (fgets(buf, sizeof(buf), input)) {
		// fprintf(stderr, "parsing: %s\n", buf);
		if (buf[0] == '\n')
			break;

		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			P.numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			P.T = dval;
		if (sscanf(buf, "\"roadwidth\":%lf", &dval) == 1)
			P.roadwidth = dval;
		if (sscanf(buf, "\"safety1\":%lf", &dval) == 1)
			P.safety1 = dval;
		if (sscanf(buf, "\"safety2\":%lf", &dval) == 1)
			P.safety2 = dval;
		if (sscanf(buf, "\"latmar\":%lf", &dval) == 1)
			P.latmar = dval;
		if (sscanf(buf, "\"penalty(%d)\":%lf", &ival, &dval) == 2)
			penalty[ival] = dval;

		/* read initial states */
		if (sscanf(buf, "\"vdx(%d)\":%lf", &ival, &dval) == 2)
			P.vdx[ival] = dval;
		if (sscanf(buf, "\"vdy(%d)\":%lf", &ival, &dval) == 2)
			P.vdy[ival] = dval;
		if (sscanf(buf, "\"x(%d)\":%lf", &ival, &dval) == 2) {
			P.x0[ival] = dval;
			P.x_len = MAX(P.x_len, ival + 1);
		}
		if (sscanf(buf, "\"y(%d)\":%lf", &ival, &dval) == 2) {
			P.y0[ival] = dval;
			P.y_len = MAX(P.y_len, ival + 1);
		}
		if (sscanf(buf, "\"vx(%d)\":%lf", &ival, &dval) == 2) {
			P.vx0[ival] = dval;
			P.vx_len = MAX(P.vx_len, ival + 1);
		}
		if (sscanf(buf, "\"vy(%d)\":%lf", &ival, &dval) == 2) {
			P.vy0[ival] = dval;
			P.vy_len = MAX(P.vy_len, ival + 1);
		}
		if (sscanf(buf, "\"ax(%d)\":%lf", &ival, &dval) == 2) {
			P.ax0[ival] = dval;
			P.ax_len = MAX(P.ax_len, ival + 1);
		}
		if (sscanf(buf, "\"ay(%d)\":%lf", &ival, &dval) == 2) {
			P.ay0[ival] = dval;
			P.ay_len = MAX(P.ay_len, ival + 1);
		}

		if (sscanf(buf, "\"ego_l(%d)\":%lf", &ival, &dval) == 2) {
			P.ego_l[ival] = dval;
			P.ego_l_len = MAX(P.ego_l_len, ival + 1);
		}
		if (sscanf(buf, "\"ego_w(%d)\":%lf", &ival, &dval) == 2) {
			P.ego_w[ival] = dval;
			P.ego_w_len = MAX(P.ego_w_len, ival + 1);
		}
		if (sscanf(buf, "\"leader(%d)\":%d", &ival, &ival2) == 2) {
			P.isLeader[ival] = ival2;
			if (ival2 == 1)
				P.leaderID = ival;
		}
		
		if (sscanf(buf, "\"cmux\":%lf", &dval) == 1)
			P.cmux = dval;
		if (sscanf(buf, "\"cmuy\":%lf", &dval) == 1)
			P.cmuy = dval;
		if (sscanf(buf, "\"p1\":%d", &ival) == 1)
			p1 = ival;
		if (sscanf(buf, "\"p2\":%d", &ival) == 1)
			p2 = ival;
		if (sscanf(buf, "\"p3\":%d", &ival) == 1)
			p3 = ival;

		if (sscanf(buf, "\"umin1\":%lf", &dval) == 1)
			P.ux_min = dval;
		if (sscanf(buf, "\"umin2\":%lf", &dval) == 1)
			P.uy_min = dval;
		if (sscanf(buf, "\"umax1\":%lf", &dval) == 1)
			P.ux_max = dval;
		if (sscanf(buf, "\"umax2\":%lf", &dval) == 1)
			P.uy_max = dval;


		if (sscanf(buf, "\"obst_x(%d)\":%lf", &ival, &dval) == 2) {
			P.obst_x[ival] = dval;
			P.obst_x_len = MAX(P.obst_x_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_y(%d)\":%lf", &ival, &dval) == 2) {
			P.obst_y[ival] = dval;
			P.obst_y_len = MAX(P.obst_y_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_vx(%d)\":%lf", &ival, &dval) == 2) {
			P.obst_vx[ival] = dval;
			P.obst_vx_len = MAX(P.obst_vx_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_vy(%d)\":%lf", &ival, &dval) == 2) {
			P.obst_vy[ival] = dval;
			P.obst_vy_len = MAX(P.obst_vy_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_l(%d)\":%lf", &ival, &dval) == 2){
			P.obst_l[ival] = dval;
			P.obst_l_len = MAX(P.obst_l_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_w(%d)\":%lf", &ival, &dval) == 2){
			P.obst_w[ival] = dval;
			P.obst_w_len = MAX(P.obst_w_len, ival + 1);
		}

		/* read FDA parameters */
		if (sscanf(buf, "\"FDA_iter\":%d", &ival) == 1)
			AP.iter = ival;
		if (sscanf(buf, "\"FDA_accur\":%lf", &dval) == 1)
			AP.accur = dval;
		if (sscanf(buf, "\"FDA_a\":%lf", &dval) == 1)
			AP.a = dval;
		if (sscanf(buf, "\"FDA_sigma\":%lf", &dval) == 1)
			AP.sigma = dval;
		if (sscanf(buf, "\"FDA_rho\":%lf", &dval) == 1)
			AP.rho = dval;
		if (sscanf(buf, "\"FDA_sigma2\":%lf", &dval) == 1)
			AP.sigma2 = dval;
		if (sscanf(buf, "\"FDA_tau\":%lf", &dval) == 1)
			AP.tau = dval;
		if (sscanf(buf, "\"FDA_lin_it\":%d", &ival) == 1)
			AP.lin_it = ival;
		if (sscanf(buf, "\"FDA_restart_it\":%d", &ival) == 1)
			AP.restart_it = ival;
		
		/* rprop */
		if (sscanf(buf, "\"FDA_eta_minus\":%lf", &dval) == 1)
			AP.eta_minus = dval;
		if (sscanf(buf, "\"FDA_eta_plus\":%lf", &dval) == 1)
			AP.eta_plus = dval;
		if (sscanf(buf, "\"FDA_delta_0(%d)\":%lf", &ival, &dval) == 2)
			AP.delta_u_0[ival] = dval;
		if (sscanf(buf, "\"FDA_delta_min\":%lf", &dval) == 1)
			AP.delta_u_min = dval;
		if (sscanf(buf, "\"FDA_delta_max\":%lf", &dval) == 1)
			AP.delta_u_max = dval;

		memset(buf, 0, sizeof(buf));
	}
	if (!(P.x_len == P.y_len && P.x_len == P.vx_len && P.x_len == P.vy_len && P.ego_l_len == P.ego_w_len && P.ax_len == P.ay_len)) {
		fputs("incomplete input\n", stderr);
		exit(1);
	}
	P.n = P.x_len;

	if (!(P.obst_x_len == P.obst_y_len && P.obst_x_len == P.obst_vx_len && P.obst_x_len == P.obst_vy_len && P.obst_l_len == P.obst_w_len)) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	P.n_obs = P.obst_x_len;
}

static void printsolSim(struct vehicle *veh, int numsteps, double T) {
	
	FILE* fd2 = NULL;
	if ((fd2 = fopen("sim/data/sim.js", "w")) == NULL) {
		// printf("\nCouldn't read data file \n");
	}
	double max_pos = 0.0;
	fputs("sim = {\n", fd2);
	fputs("x: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++) {
			fprintf(fd2, "%.4f,", veh[idVeh].states[0][k]);
			if (veh[idVeh].states[0][k] > max_pos)
				max_pos = veh[idVeh].states[0][k];
		}
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", obs.x[idVeh][0][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("ax: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].control[0][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("ay: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].control[1][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("y: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].states[1][k]);
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", obs.x[idVeh][1][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vx: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].states[2][k]);
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", obs.v[idVeh][0][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vdx: [", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fd2, "%.4f,", P.vdx[idVeh]);
	fputs("],\n", fd2);

	fputs("id: [", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fd2, "%d,", idVeh + 1);
	fputs("],\n", fd2);

	fputs("Cx: [", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fd2, "%.4f,", P.ego_l[idVeh]);
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++)
		fprintf(fd2, "%.4f,", P.obst_l[idVeh]);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fd2, "%.4f,", P.ego_w[idVeh] );
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++)
		fprintf(fd2, "%.4f,", P.obst_w[idVeh]);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", P.n + P.n_obs);
	fprintf(fd2, "k:%d,\n", P.numsteps);
	fprintf(fd2, "roadlength:%f,\n", max_pos);
	fprintf(fd2, "roadwidth:%f,\n", P.roadwidth);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", P.T);
	fprintf(fd2, "safety1:%f,\n", P.safety1);
	fprintf(fd2, "safety2:%f,\n", P.safety2);
	fprintf(fd2, "roadbound:%f,\n", 1.2);

	fprintf(fd2, "attraction_x:%f,\n", at.len);
	fprintf(fd2, "attraction_y:%f,\n", at.width);
	fprintf(fd2, "leaderID:%d,\n", P.leaderID);

	fprintf(fd2, "offset_count:%d,\n", 0);
	fprintf(fd2, "limits: [],\n");
	fprintf(fd2, "slopes: [],\n");
	fprintf(fd2, "offsets: [],\n");
	fprintf(fd2, "};\n\n");

	fclose(fd2);
}

static void plotFile(struct vehicle *veh, int numsteps, double T) {
	
	FILE *fout;
	fout = fopen("output.py", "w");

	/* python file */
	fputs("sim = {\n", fout);
	fputs("'x': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].states[0][k]);
		fputs("],\n", fout);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", obs.x[idVeh][0][k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	
	double uxTrans;
	// fputs("'ax': [\n", fout);
	// for (int idVeh = 0; idVeh < P.n; idVeh++) {
	// 	fputs("\t[", fout);
	// 	for (int k = 0; k < P.numsteps; k++)
	// 		fprintf(fout, "%.4f,", veh[idVeh].control[0][k]);
	// 	fputs("],\n", fout);
	// }
	// fputs("],\n", fout);

	fputs("'ax': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++) {
			uxTrans = (1.0 - veh[idVeh].control[0][k]) * veh[idVeh].uxLB[k] + veh[idVeh].control[0][k] * veh[idVeh].uxUB[k];
			fprintf(fout, "%.4f,", uxTrans);
		}
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'axLB': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].uxLB[k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'axUB': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].uxUB[k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	double uyTrans;
	// fputs("'ay': [\n", fout);
	// for (int idVeh = 0; idVeh < P.n; idVeh++) {
	// 	fputs("\t[", fout);
	// 	for (int k = 0; k < P.numsteps; k++)
	// 		fprintf(fout, "%.4f,", veh[idVeh].control[1][k]);
	// 	fputs("],\n", fout);
	// }
	// fputs("],\n", fout);

	fputs("'ay': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++) {
			uyTrans = (1.0 - veh[idVeh].control[1][k]) * veh[idVeh].uyLB[k] + veh[idVeh].control[1][k] * veh[idVeh].uyUB[k];
			fprintf(fout, "%.4f,", uyTrans);
		}
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'ayLB': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].uyLB[k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'ayUB': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].uyUB[k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'y': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].states[1][k]);
		fputs("],\n", fout);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", obs.x[idVeh][1][k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'vx': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].states[2][k]);
		fputs("],\n", fout);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", obs.v[idVeh][0][k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'vy': [\n", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", veh[idVeh].states[3][k]);
		fputs("],\n", fout);
	}
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++) {
		fputs("\t[", fout);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fout, "%.4f,", obs.v[idVeh][1][k]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'vdx': [", fout);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fout, "%.4f,", P.vdx[idVeh]);
	fputs("],\n", fout);

	fprintf(fout, "'n':%d,\n", P.n + P.n_obs);
	fprintf(fout, "'n_obs':%d,\n", P.n_obs);
	fprintf(fout, "'k':%d,\n", P.numsteps);
	fprintf(fout, "'Step':%f,\n", P.T);
	fprintf(fout, "}\n\n");


	fclose(fout);
}

int main() {

	read_P();

	P.ux_min = -2.0;
	P.ux_max = 1.5;


	at.len = 50.0;
	at.width = 1.0;
	at.p1 = 24; at.p2 = 18;

	clock_t start, start_FDA, end, end_FDA;
	double cpu_time_used = 0, cpu_time_used_FDA = 0;
	start = clock();
	
	allocations(P.numsteps);
	obst_trajectories(P.numsteps, P.T);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		veh[idVeh].states[0][0] = P.x0[idVeh];
		veh[idVeh].states[1][0] = P.y0[idVeh];
		veh[idVeh].states[2][0] = P.vx0[idVeh];
		veh[idVeh].states[3][0] = P.vy0[idVeh];
	}

	if (SDIR == 5) {
		for (int idVeh = 0; idVeh < P.n; idVeh++) {
			for (int i = 0; i < NU; i++) {
				for (int j = 0; j < P.numsteps; j++) {
					veh[idVeh].delta_u[i][j] = AP.delta_u_0[i];
				}
			}
		}
	}

	double lm = 0, yLB, yUB, vyLB, vyUB;
	double vx0, y0, vy0;
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		lm = P.ego_w[idVeh] / 2.0;
		vx0 = P.vx0[idVeh];
		y0 = P.y0[idVeh];
		vy0 = P.vy0[idVeh];
		
		veh[idVeh].yLim[0][0] = lm;
		veh[idVeh].yLim[1][0] = P.roadwidth - lm;
		veh[idVeh].vyLim[0][0] = veh[idVeh].vyLim[1][0] = 0.0;

		for (int i = 0; i < P.numsteps; i++) {
			veh[idVeh].vyLim[0][i + 1] = vyLB = veh[idVeh].vyLim[0][0];
			veh[idVeh].vyLim[1][i + 1] = vyUB = veh[idVeh].vyLim[1][0];
			veh[idVeh].yLim[0][i + 1]  = yLB  = veh[idVeh].yLim[0][0];
			veh[idVeh].yLim[1][i + 1]  = yUB  = veh[idVeh].yLim[1][0];
			
			veh[idVeh].control[0][i] = (-LONK3 * P.vx0[idVeh]) / (-P.ux_max - LONK3 * P.vx0[idVeh]);
			veh[idVeh].control[1][i] = ((LATK1 * (-P.y0[idVeh] + yLB) + LATK2 * (vyLB - P.vy0[idVeh]))
				/ ((LATK1 * (-P.y0[idVeh] + yLB) + LATK2 * (vyLB - P.vy0[idVeh])) - (LATK1 * (-P.y0[idVeh] + yUB) + LATK2 * (vyUB - P.vy0[idVeh]))) );
		}
	}

	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_states(veh, veh[idVeh].states, veh[idVeh].control, P.numsteps, P.T, veh[idVeh].states, idVeh); }

	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		for (int i = 0; i < NX; i++) {
			for (int k = 0; k < P.numsteps + 1; k++) {
				states_all_veh[idVeh][i][k] = veh[idVeh].states[i][k];
			}
		}
	}
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		for (int i = 0; i < NU; i++) {
			for (int k = 0; k < P.numsteps; k++) {
				control_all_veh[idVeh][i][k] = veh[idVeh].control[i][k];
			}
		}
	}
	
	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_costates(veh, states_all_veh, control_all_veh, P.numsteps, P.T, veh[idVeh].costates, idVeh); }
	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_reduced_gradient(veh, states_all_veh, control_all_veh, P.numsteps, P.T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh); }

	start_FDA = clock();
	optimization(veh, criterion, P.numsteps, P.T, AP.iter, AP.accur, AP.a, AP.sigma, AP.rho, AP.tau, AP.lin_it, AP.restart_it);
	end_FDA = clock();
	cpu_time_used_FDA = ((double)(end_FDA - start_FDA)) / CLOCKS_PER_SEC;
	fprintf(stderr, "CPU Time: %.4f\n", cpu_time_used_FDA);
	
	printsolSim(veh, P.numsteps, P.T);
	plotFile(veh, P.numsteps, P.T);

	return 0;
}