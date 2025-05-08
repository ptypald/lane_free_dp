#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>

#define MAX(a,b) (a > b)?(a):(b)
#define MIN(a,b) (a < b)?(a):(b)
#define MYSIGN(a) (a > 0.0)?(1.0):(-1.0)
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
#define EPSILONW 0.1

double LATK1 = 0.25 / 1;
double LATK2 = 2.0 * sqrt(LATK1) - LATK1 * (0.25 / 2.0);

double LONK3 = 4.0;
double BETA = 0.03;

//                        { UX,  UY, VDX, VDY, COL-AVD1, DIS-U, VCP, COL-AVD2, FLBL/FLBU, DIS-L };
static double penalty[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

double c_min[] = { 0.0, 0.0 };
double c_max[] = { 1.0, 1.0 };

/* ellipse parameters */
int p1 = 2, p2 = 2, p3 = 2;

static double* grad_norm;
static double*** states_all_veh, *** control_all_veh;
/* debugging */
int ls_no_iter = 0, brac_it = 0, sect_it = 0;

struct vehicle {
	double** states, ** costates;
	double** control;
	double** reduced_g;
	double** reduced_g_prev;
	double** s_dir, ** temp_states, ** temp_costates, ** temp_control, ** temp_reduced_g;
	double** control_prev, ** control_cur, ** reduced_g_cur;
	double** control_a1, ** control_a2;
	double** reduced_g_a1, ** reduced_g_a2;
	double** control_alpha, ** reduced_g_alpha;
	double** states_prev, ** states_cur, ** costates_cur;
	double** states_a1, ** costates_a1;
	double** states_a2, ** costates_a2;
	double** states_alpha, ** costates_alpha;
	double** grad_dif;

	// double **x1Limits, **x1Limitsv, **x2Limits, **x2Limitsv, **x2Limitsa;
	double** xLim, ** vxLim, ** yLim, ** vyLim;
	double* uxLB, * uxUB, * uyLB, * uyUB;

	/* rprop */
	double** delta_u;
	double** rprop_best_control;
};

struct vehicle* veh;

static const char* criterion = "gradient";
static struct {
	int iter;
	double accur;
	double a;
	double sigma, rho, sigma2, tau;
	int lin_it, restart_it;

	/* rprop */
	double eta_minus, eta_plus;
	double delta_u_0;
	double delta_u_min, delta_u_max;
}AP = { 0 };

static struct {
	double* vdx, * vdy, T, * ego_l, * ego_w;
	double* x0, * y0, * vx0, * vy0, * ax0, * ay0;

	int x_len, y_len, vx_len, vy_len, ego_l_len, ego_w_len, ax_len, ay_len;

	int numsteps, n, rad;
	double safety, roadwidth, latmar, elc, slp;

	double cmux, cmuy;
	double ux_min, ux_max, uy_min, uy_max;

	double* obst_x, * obst_y, * obst_vx, * obst_vy, * obst_l, * obst_w;
	int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len, obst_l_len, obst_w_len, obst_ax_len, obst_ay_len;
	int n_obs;
}P = { 0 };

static struct {
	bool apply;
}r = { 0 };

static struct {
	double*** x;
	double*** v;
}obs = { 0 };

static double inner_product(double** A, double** B, int numsteps, int nu) {
	int i, j;
	double result;
	result = 0.0;

	for (i = 0; i < nu; i++) {
		for (j = 0; j < numsteps; j++) {
			result += A[i][j] * B[i][j];
		}
	}
	return result;
}

static void inner_product_1(double** A, double** B, int numsteps, int nu, double* result) {

	for (int i = 0; i < nu; i++) {
		for (int j = 0; j < numsteps; j++) {
			result[i] += A[i][j] * B[i][j];
		}
	}
}

static double sign(double a) {
	if (a > 0.0) {
		return 1.0;
	}
	else if (a < 0.0) {
		return -1.0;
	}
	else {
		return 0.0;
	}
}

static bool restart(struct vehicle* veh, int numsteps, int nu, int iter) {
	int idVeh;
	double inner_1 = 0.0, inner_2 = 0.0, inner_3 = 0.0;

	for (idVeh = 0; idVeh < P.n; idVeh++) {
		inner_1 += inner_product(veh[idVeh].s_dir, veh[idVeh].reduced_g, numsteps, nu);
		inner_2 += inner_product(veh[idVeh].s_dir, veh[idVeh].s_dir, numsteps, nu);
		inner_3 += inner_product(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, nu);
	}

	if (inner_1 >= -AP.sigma2 * sqrt(inner_2 * inner_3)) {
		fprintf(stderr, "--restart true: %.7f >= %.7f\n", inner_1, -AP.sigma2 * sqrt(inner_2 * inner_3));
		return true;
	}

	return false;
}

static double quadratic_interpolation(double a, double b, double dfa, double dfb) {
	double alpha;

	alpha = (b * dfa - a * dfb) / (dfa - dfb);

	return alpha;
}

static double cubic_interpolation(double a, double b, double fa, double fb, double dfa, double dfb) {
	double z, w;
	double alpha;

	z = 3.0 * ((fa - fb) / (b - a)) + dfa + dfb;
	w = pow(z, 2) - dfa * dfb;

	if (w > 0) {
		w = sqrt(w);
		alpha = b - (((dfb + w - z) / (dfb - dfa + 2.0 * w)) * (b - a));
	}
	else {
		//fprintf(stderr, "negative sqrt in cubic interpolation, apply bisection \n");
		alpha = (a + b) / 2.0;
	}
	return alpha;
}

static bool termination(struct vehicle* veh, const char* termination, double accuracy, double cost_cur, double cost_prev, int numsteps) {
	double temp = 0.0;

	if (strncmp(termination, "relative", strlen(termination)) == 0) {
		temp = (cost_cur - cost_prev) / (cost_cur + 0.001);
		fprintf(stderr, "%f, %f, %f, %f\n", cost_cur, cost_prev, temp, accuracy);
		if (fabs(temp) < accuracy) {
			fprintf(stderr, "returned true\n");
			return true;
		}
	}
	else if (strncmp(termination, "gradient", strlen(termination)) == 0) {
		for (int idVeh = 0; idVeh < P.n; idVeh++)
			temp += inner_product(veh[idVeh].temp_reduced_g, veh[idVeh].temp_reduced_g, numsteps, NU);
		if (temp < accuracy) {
			return true;
		}
	}
	return false;
}

static bool termination_1(struct vehicle* veh, const char* termination, double accuracy, double cost_cur, double cost_prev, int numsteps, double* grad_norm, int it) {
	double temp;

	if (strncmp(termination, "relative", strlen(termination)) == 0) {
		temp = (cost_cur - cost_prev) / (cost_cur + 0.001);
		if (fabs(cost_cur - cost_prev) < 1e-12) {
			return true;
		}
	}
	else if (strncmp(termination, "gradient", strlen(termination)) == 0) {
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
			//obs.x[i][1][j + 1] = obs.x[i][1][j];

			//this is working-sinusoidal-y-trajectory-synchronous
			//obs.x[i][1][j + 1] = 9 * sin(0.05 * j) + 10;  // roadwidth-20
			//obs.x[i][1][j + 1] = 6.5 * sin(0.08 * j) + 7.5; // roadwidth-15

			//has crashes-sinusoidal-y-trajectory
			//obs.x[i][1][j + 1] = 6.1 * sin(0.3 * j) + 7.1;

			//normal-y-trajectory
			obs.x[i][1][j + 1] = obs.x[i][1][j];

			//asynchronous-obstacle-traj
			//obs.x[i][1][j + 1] = 6.1 * sin(0.04 * j + i) + 7.1;
		}
		// }
	}
}

static void compute_states(struct vehicle* veh, double** x_init, double*** states, double** control, int numsteps, double T, double** states_calc, int idVeh) {
	int i;

	for (i = 0; i < NX; i++) {
		states_calc[i][0] = x_init[i][0];
	}

	double epsilon = 0.1, UXMIN = P.ux_min;
	double uxTrans, uyTrans, yLB, yUB, vyLB, vyUB;
	for (i = 0; i < numsteps; i++) {

		yLB = veh[idVeh].yLim[0][i]; yUB = veh[idVeh].yLim[1][i]; vyLB = veh[idVeh].vyLim[0][i]; vyUB = veh[idVeh].vyLim[1][i];

		veh[idVeh].uxLB[i] = 0.5 * (UXMIN - LONK3 * states_calc[2][i] + Sqrt(epsilon + Power(-UXMIN - LONK3 * states_calc[2][i], 2)));
		veh[idVeh].uxUB[i] = P.ux_max;

		veh[idVeh].uyLB[i] = -LATK1 * (states_calc[1][i] - yLB) - LATK2 * (states_calc[3][i] - vyLB);
		veh[idVeh].uyUB[i] = -LATK1 * (states_calc[1][i] - yUB) - LATK2 * (states_calc[3][i] - vyUB);


		//fprintf(stderr, "(%d, %d) %.4f \t %.4f \t %.4f \t %.4f \n", idVeh, i, veh[idVeh].uxLB[i], veh[idVeh].uxUB[i], veh[idVeh].uyLB[i], veh[idVeh].uyUB[i]);


		uxTrans = (1.0 - control[0][i]) * veh[idVeh].uxLB[i] + control[0][i] * veh[idVeh].uxUB[i];
		uyTrans = (1.0 - control[1][i]) * veh[idVeh].uyLB[i] + control[1][i] * veh[idVeh].uyUB[i];

		states_calc[0][i + 1] = states_calc[0][i] + states_calc[2][i] * T + 0.5 * uxTrans * pow(T, 2);
		states_calc[1][i + 1] = states_calc[1][i] + states_calc[3][i] * T + 0.5 * uyTrans * pow(T, 2);
		states_calc[2][i + 1] = states_calc[2][i] + uxTrans * T;
		states_calc[3][i + 1] = states_calc[3][i] + uyTrans * T;
	}
}

static void compute_costates(struct vehicle* veh, double*** states, double*** control, int numsteps, double T, double** costates_calc, double** reduced_g_calc, int idVeh) {
	double dvdx, dvdy;
	double x, y, vx, vy, ux, uy;

	double duxdx, duxdvx, duydy, duydvy;

	double l_m = P.ego_w[idVeh] / 2.0;
	double y0 = P.roadwidth - l_m;

	double epsilon = 0.1, UXMAX = P.ux_max, UXMIN = P.ux_min;

	double p = 6.0, gamma = 1.0;

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

		double costate_term_1 = 0.0, costate_term_2 = 0.0, costate_term_3 = 0.0;
		double L = 5.0, w = 2.0, tau = P.safety;

		/* collision avoidance term derivatice w.r.t x, y, v */
		for (int ob = 0; ob < P.n_obs; ob++) {
			double obsvx = obs.v[ob][0][i], obsvy = obs.v[ob][1][i], obsx = obs.x[ob][0][i], obsy = obs.x[ob][1][i];
			double x = states[idVeh][0][i], y = states[idVeh][1][i], vx = states[idVeh][2][i], vy = states[idVeh][3][i];
			
			//Elliptical Collision Avoidance
			if (x <= obsx) {
				costate_term_1 += (-2*(-obsx + x))/(Power(L + tau*vx,2)*
     				Power(1 + Power(-obsx + x,2)/Power(L + tau*vx,2) + Power(-obsy + y,2)/Power(w,2),2));
				costate_term_2 += (-2*(-obsy + y))/(Power(w,2)*Power(1 + Power(-obsx + x,2)/Power(L + tau*vx,2) + 
       				Power(-obsy + y,2)/Power(w,2),2));
				costate_term_3 += (2*tau*Power(-obsx + x,2))/
   					(Power(L + tau*vx,3)*Power(1 + Power(-obsx + x,2)/Power(L + tau*vx,2) + 
       				Power(-obsy + y,2)/Power(w,2),2));
			}
			else {
				costate_term_1 += (-2*(-obsx + x))/(Power(L + obsvx*tau,2)*
     				Power(1 + Power(-obsx + x,2)/Power(L + obsvx*tau,2) + Power(-obsy + y,2)/Power(w,2),2));
				costate_term_2 += (-2*(-obsy + y))/(Power(w,2)*Power(1 + Power(-obsx + x,2)/Power(L + tau*vx,2) + 
       				Power(-obsy + y,2)/Power(w,2),2));
				costate_term_3 += 0.0;
			}
		}

		dvdx = penalty[2] * (states[idVeh][2][i] - P.vdx[idVeh]) / pow(P.vdx[idVeh], 2);
		dvdy = penalty[3] * (states[idVeh][3][i] - P.vdy[idVeh]) / pow(2.0, 2);
		
		double dudx, dudvx, dudy, dudvy;
		double yLB = veh[idVeh].yLim[0][i], yUB = veh[idVeh].yLim[1][i], vyLB = veh[idVeh].vyLim[0][i], vyUB = veh[idVeh].vyLim[1][i];

		dudx = 0.0;
		dudy = 1. * penalty[1] * (-(LATK1 * (1 - uy)) - LATK1 * uy) * ((1 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) + uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB))) / pow(P.uy_max, 2);
		dudvx = 0.5 * penalty[0] * (1 - ux) * (-LONK3 - (LONK3 * (-UXMIN - LONK3 * vx)) / Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2))) * (ux * UXMAX + 0.5 * (1 - ux) * (UXMIN - LONK3 * vx + Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2)))) / pow(P.ux_max, 2);
		dudvy = 1. * penalty[1] * (-(LATK2 * (1 - uy)) - LATK2 * uy) * ((1 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) + uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB))) / pow(P.uy_max, 2);

		djxdvx = ...;
		djydy = ...;
		djydvy = ...;


		/* costates calculation */
		costates_calc[0][i] = costates_calc[0][i + 1] + penalty[4] * costate_term_1;
		costates_calc[1][i] = djydy + penalty[4] * costate_term_2 + dudy + costates_calc[3][i + 1] * T * (-(LATK1 * (1 - uy)) - LATK1 * uy) + costates_calc[1][i + 1] * (1 + 0.5 * Power(T, 2) * (-(LATK1 * (1 - uy)) - LATK1 * uy));
		costates_calc[2][i] = penalty[4] * costate_term_3 + dvdx + dudvx + costates_calc[2][i + 1] * (1 + 0.5 * T * (1 - ux) * (-LONK3 - (LONK3 * (-UXMIN - LONK3 * vx)) / Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2)))) + costates_calc[0][i + 1] * (T + 0.25 * Power(T, 2) * (1 - ux) * (-LONK3 - (LONK3 * (-UXMIN - LONK3 * vx)) / Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2))));
		costates_calc[3][i] = dvdy + dudvy + costates_calc[3][i + 1] * (1 + T * (-(LATK2 * (1 - uy)) - LATK2 * uy)) + costates_calc[1][i + 1] * (T + 0.5 * Power(T, 2) * (-(LATK2 * (1 - uy)) - LATK2 * uy));

		i = i - 1;
	}
}

static double compute_cost(struct vehicle* veh, double*** states, double*** control, int numsteps, double T, int idVeh) {
	int i, j;
	double ux, uy, vdx, vdy, ellipse_cost, vcp;
	
	double uxTrans, uyTrans;
	double cost = 0.0;

	double l_m = P.ego_w[idVeh] / 2.0;
	double y0 = P.roadwidth - l_m;

	double UXMAX = P.ux_max, UXMIN = P.ux_min, epsilon = 0.1;
	double yLB, yUB, vyLB, vyUB;

	double p = 6.0, gamma = 1.0;
	double L = 5.0, w = 2.0, tau = P.safety;
	
	for (i = 0; i < numsteps; i++) {
		ellipse_cost = 0.0;

		/* collision avoidance term */
		for (int ob = 0; ob < P.n_obs; ob++) {
			double obsvx = obs.v[ob][0][i], obsvy = obs.v[ob][1][i], obsx = obs.x[ob][0][i], obsy = obs.x[ob][1][i];
			double x = states[idVeh][0][i], y = states[idVeh][1][i], vx = states[idVeh][2][i], vy = states[idVeh][3][i];
			 
			//Eliptical-distance Collision Avoidance
			if (x <= obsx)
				ellipse_cost += 1/(1 + Power(-obsx + x,2)/Power(L + tau*vx,2) + Power(-obsy + y,2)/Power(w,2));
			else
				ellipse_cost += 1/(1 + Power(-obsx + x,2)/Power(L + obsvx*tau,2) + Power(-obsy + y,2)/Power(w,2));

		}
		ellipse_cost *= penalty[4];

		yLB = veh[idVeh].yLim[0][i]; yUB = veh[idVeh].yLim[1][i]; vyLB = veh[idVeh].vyLim[0][i]; vyUB = veh[idVeh].vyLim[1][i];

		uxTrans = control[idVeh][0][i] * UXMAX + 0.5 * (1.0 - control[idVeh][0][i]) * (UXMIN - LONK3 * states[idVeh][2][i] + Sqrt(epsilon + Power(-UXMIN - LONK3 * states[idVeh][2][i], 2)));
		uyTrans = (1.0 - control[idVeh][1][i]) * (-(LATK2 * (states[idVeh][3][i] - vyLB)) - LATK1 * (states[idVeh][1][i] - yLB)) + control[idVeh][1][i] * (-(LATK2 * (states[idVeh][3][i] - vyUB)) - LATK1 * (states[idVeh][1][i] - yUB));

		/* terms for u (longtitudinal/lateral) and desired speed */
		ux = 0.5 * penalty[0] * pow(uxTrans, 2) / pow(P.ux_max, 2);
		uy = 0.5 * penalty[1] * pow(uyTrans, 2) / pow(P.uy_max, 2);
		vdx = 0.5 * penalty[2] * pow((states[idVeh][2][i] - P.vdx[idVeh]), 2) / pow(P.vdx[idVeh], 2);
		vdy = 0.5 * penalty[3] * pow((states[idVeh][3][i] - P.vdy[idVeh]), 2) / pow(2.0, 2);

		cost += ux + uy + vdx + vdy + ellipse_cost;
	}
	
	return cost;
}

static void compute_reduced_gradient(struct vehicle* veh, double*** states, double*** control, int numsteps, double T, double** costates, double** reduced_g_calc, int idVeh) {
	int i, j;
	double x, y, vx, vy, ux, uy;
	double epsilon = 0.1, UXMAX = P.ux_max, UXMIN = P.ux_min;

	for (i = 0; i < numsteps; i++) {

		x = states[idVeh][0][i]; y = states[idVeh][1][i];
		vx = states[idVeh][2][i]; vy = states[idVeh][3][i];
		ux = control[idVeh][0][i]; uy = control[idVeh][1][i];

		double yLB = veh[idVeh].yLim[0][i], yUB = veh[idVeh].yLim[1][i], vyLB = veh[idVeh].vyLim[0][i], vyUB = veh[idVeh].vyLim[1][i];

		reduced_g_calc[0][i] = penalty[0] * (UXMAX - 0.5 * (UXMIN - LONK3 * vx + Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2)))) * (ux * UXMAX + 0.5 * (1 - ux) * (UXMIN - LONK3 * vx + Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2)))) / pow(P.ux_max, 2) + 
			costates[2][i + 1] * T * (UXMAX - 0.5 * (UXMIN - LONK3 * vx + Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2)))) + 
			0.5 * costates[0][i + 1] * Power(T, 2) * (UXMAX - 0.5 * (UXMIN - LONK3 * vx + Sqrt(epsilon + Power(-UXMIN - LONK3 * vx, 2))));
		reduced_g_calc[1][i] = penalty[1] * ((1 - uy) * (-(LATK2 * (vy - vyLB)) - LATK1 * (y - yLB)) + uy * (-(LATK2 * (vy - vyUB)) - LATK1 * (y - yUB))) * (LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB)) / pow(P.uy_max, 2) + 
			costates[3][i + 1] * T * (LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB)) + 
			0.5 * costates[1][i + 1] * Power(T, 2) * (LATK2 * (vy - vyLB) - LATK2 * (vy - vyUB) + LATK1 * (y - yLB) - LATK1 * (y - yUB));
	}
	for (i = 0; i < NU; i++) {
		for (j = 0; j < numsteps; j++) {
			if ((control[idVeh][i][j] == c_min[i]) && (reduced_g_calc[i][j] > 0.0)) {
				reduced_g_calc[i][j] = 0.0;
			}
			else if ((control[idVeh][i][j] == c_max[i]) && (reduced_g_calc[i][j] < 0.0)) {
				reduced_g_calc[i][j] = 0.0;
			}
		}
	}

}

static void control_bound(double** control, double* control_min, double* control_max, int numsteps, int nu, int idVeh) {
	int i, j;

	for (i = 0; i < nu; i++) {
		for (j = 0; j < numsteps; j++) {
			if (control[i][j] < control_min[i]) {
				control[i][j] = control_min[i];
			}
			else if (control[i][j] > control_max[i]) {
				control[i][j] = control_max[i];
			}
		}
	}
}

static double line_search(struct vehicle* veh, double a, double sigma, double tau, double rho, int lin_it, int numsteps, double T, int iter, double cost_prev_it) {
	int i, j, b_it, idVeh, k;
	double f0, df0, fai_1, fai, dfai;
	double a_cur, a_prev, a_temp, a1 = 0.0, a2 = 0.0;
	double fa1, dfa1, fa2, dfa2;
	double fa, dfa, alpha;

	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < NX; i++) {
			for (k = 0; k < P.numsteps + 1; k++) {
				states_all_veh[idVeh][i][k] = veh[idVeh].states[i][k];
			}
		}
	}
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < NU; i++) {
			for (k = 0; k < P.numsteps; k++) {
				control_all_veh[idVeh][i][k] = veh[idVeh].control[i][k];
			}
		}
	}

	f0 = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++)
		f0 += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh);
		//compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh);
	}
	for (idVeh = 0; idVeh < P.n; idVeh++)
		compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh);

	df0 = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++)
		df0 += inner_product(veh[idVeh].reduced_g, veh[idVeh].s_dir, numsteps, NU);

	/* check */
	a_cur = (-2.0 * (cost_prev_it - f0)) / df0;
	if (iter == 0 || a_cur < 0.0)
		a_cur = a;		/* a_cur: initial step, given from the user */
	a_prev = 0.0;	/* Initialization of a_prev */

	/* Bracketing */
	/*            */
	brac_it = sect_it = 0;
	for (b_it = 0; b_it < lin_it; b_it++) {
		/* debugging */
		brac_it = b_it + 1;

		/* Calculation of previous control*/
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (j = 0; j < numsteps; j++) {
					veh[idVeh].control_prev[i][j] = veh[idVeh].control[i][j] + a_prev * veh[idVeh].s_dir[i][j];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].control_prev, c_min, c_max, numsteps, NU, idVeh);

		/* Calculation of previous states */
		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control_prev, numsteps, T, veh[idVeh].states_prev, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].states_prev[i][k];
				}
			}
		}
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].control_prev[i][k];
				}
			}
		}

		/* Calculation of F(a_i-1) */
		fai_1 = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			fai_1 += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		/* Calculation of current control */
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (j = 0; j < numsteps; j++) {
					veh[idVeh].control_cur[i][j] = veh[idVeh].control[i][j] + a_cur * veh[idVeh].s_dir[i][j];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].control_cur, c_min, c_max, numsteps, NU, idVeh);

		/* Calculation of current states */
		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control_cur, numsteps, T, veh[idVeh].states_cur, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].states_cur[i][k];
				}
			}
		}
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].control_cur[i][k];
				}
			}
		}

		/* Calculation of F(a_i) */
		fai = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			fai += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		/* If one of the below conditions is satisfied, bracketing phase ends */
		if (((fai > f0 + a_cur * rho * df0) || (fai >= fai_1))) {
			a1 = a_prev;
			a2 = a_cur;

			break;
		}
		else {
			/* Otherwise */
			/* Compute current costates */
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NX; i++) {
					for (k = 0; k < P.numsteps + 1; k++) {
						states_all_veh[idVeh][i][k] = veh[idVeh].states_cur[i][k];
					}
				}
			}
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (k = 0; k < P.numsteps; k++) {
						control_all_veh[idVeh][i][k] = veh[idVeh].control_cur[i][k];
					}
				}
			}

			for (idVeh = 0; idVeh < P.n; idVeh++) {
				compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_cur, veh[idVeh].reduced_g_cur, idVeh);
				//compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_cur, veh[idVeh].reduced_g_cur, idVeh);
			}
			for (idVeh = 0; idVeh < P.n; idVeh++)
				compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_cur, veh[idVeh].reduced_g_cur, idVeh);

			/* Compute F'(a_i) */
			dfai = 0.0;
			for (idVeh = 0; idVeh < P.n; idVeh++)
				dfai += inner_product(veh[idVeh].reduced_g_cur, veh[idVeh].s_dir, numsteps, NU);

			/* If the below condition is satisfied, line search phase ends */
			if (fabs(dfai) <= -sigma * df0) {
				a1 = a_cur;
				a2 = a_cur;

				// fprintf(stderr, "-- condition 2 (return at %d iters): End of line search phase with alpha: %f \n", b_it, a_cur);
				return a_cur;
			}
			else if (dfai >= 0.0) {
				a1 = a_prev;
				a2 = a_cur;

				// fprintf(stderr, "-- condition 3 (break at %d iters): End of bracketing phase with a1: %f and a2: %f \n", b_it, a1, a2);
				break;
			}
			a_temp = a_cur;
			a_cur = a_cur + tau * (a_cur - a_prev);
			a_prev = a_temp;

			a1 = a_prev;
			a2 = a_cur;
		}
	}

	/* Sectioning */
	/*            */
	int count = 0;
	while (true) {
		/* debugging */
		sect_it = count + 1;

		if (fabs(a1 - a2) < 1.e-22 || count == lin_it) {
			if (count == lin_it) {
				if (fa > f0)
					r.apply = true;
				fprintf(stderr, "Bracketing upper bound reached (f0: %.4f -- fa: %.4f)\n", f0, fa);
				return alpha;
			}
			return ((a1 + a2) / 2.0);
		}

		/* calculation for a1 */
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (j = 0; j < numsteps; j++) {
					veh[idVeh].control_a1[i][j] = veh[idVeh].control[i][j] + a1 * veh[idVeh].s_dir[i][j];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].control_a1, c_min, c_max, numsteps, NU, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control_a1, numsteps, T, veh[idVeh].states_a1, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].states_a1[i][k];
				}
			}
		}
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].control_a1[i][k];
				}
			}
		}

		fa1 = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			fa1 += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_a1, veh[idVeh].reduced_g_a1, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_a1, veh[idVeh].reduced_g_a1, idVeh);

		dfa1 = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			dfa1 += inner_product(veh[idVeh].reduced_g_a1, veh[idVeh].s_dir, numsteps, NU);

		/* calculation for a2 */
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < numsteps; k++) {
					veh[idVeh].control_a2[i][k] = veh[idVeh].control[i][k] + a2 * veh[idVeh].s_dir[i][k];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].control_a2, c_min, c_max, numsteps, NU, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control_a2, numsteps, T, veh[idVeh].states_a2, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].states_a2[i][k];
				}
			}
		}
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].control_a2[i][k];
				}
			}
		}

		fa2 = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			fa2 += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_a2, veh[idVeh].reduced_g_a2, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_a2, veh[idVeh].reduced_g_a2, idVeh);

		dfa2 = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			dfa2 += inner_product(veh[idVeh].reduced_g_a2, veh[idVeh].s_dir, numsteps, NU);

		/*  Interpolation Methods for a1 and a2 */
		alpha = cubic_interpolation(a1, a2, fa1, fa2, dfa1, dfa2);
		// alpha = quadratic_interpolation(a1, a2, dfa1, dfa2);

		/* calculation for alpha */
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (j = 0; j < numsteps; j++) {
					veh[idVeh].control_alpha[i][j] = veh[idVeh].control[i][j] + alpha * veh[idVeh].s_dir[i][j];
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].control_alpha, c_min, c_max, numsteps, NU, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control_alpha, numsteps, T, veh[idVeh].states_alpha, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NX; i++) {
				for (k = 0; k < P.numsteps + 1; k++) {
					states_all_veh[idVeh][i][k] = veh[idVeh].states_alpha[i][k];
				}
			}
		}
		for (idVeh = 0; idVeh < P.n; idVeh++) {
			for (i = 0; i < NU; i++) {
				for (k = 0; k < P.numsteps; k++) {
					control_all_veh[idVeh][i][k] = veh[idVeh].control_alpha[i][k];
				}
			}
		}

		fa = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			fa += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		// fprintf(stderr, "------ a1: %.4f \t a2: %.4f \t fa1: %.4f \t fa2: %.4f \t dfa1: %.4f \t dfa2: %.4f \t alpha: %.4f \t fa: %.4f \n", a1, a2, fa1, fa2, dfa1, dfa2, alpha, fa);
		if (fa > f0 + rho * alpha * df0 || fa >= fa1) {
			a2 = alpha;
		}
		else {
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NX; i++) {
					for (k = 0; k < P.numsteps + 1; k++) {
						states_all_veh[idVeh][i][k] = veh[idVeh].states_alpha[i][k];
					}
				}
			}
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (k = 0; k < P.numsteps; k++) {
						control_all_veh[idVeh][i][k] = veh[idVeh].control_alpha[i][k];
					}
				}
			}

			for (idVeh = 0; idVeh < P.n; idVeh++)
				compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_alpha, veh[idVeh].reduced_g_alpha, idVeh);

			for (idVeh = 0; idVeh < P.n; idVeh++)
				compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].costates_alpha, veh[idVeh].reduced_g_alpha, idVeh);

			dfa = 0.0;
			for (idVeh = 0; idVeh < P.n; idVeh++)
				dfa += inner_product(veh[idVeh].reduced_g_alpha, veh[idVeh].s_dir, numsteps, NU);

			if (fabs(dfa) <= -sigma * df0) {
				// fprintf(stderr, "-- condition (return at %d iters): End of sectioning phase with alpha: %4f \n", count, alpha);
				return alpha;
			}
			else if (dfa > 0.0) {
				a2 = alpha;
			}
			else if (dfa < 0.0) {
				a1 = alpha;
			}
		}
		count++;
	}
}

static void fletcher_reeves(struct vehicle* veh, int numsteps, int nu) {
	int i, j, idVeh;
	double beta;

	double nom = 0.0, denom = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		nom += inner_product(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, nu);
		denom += inner_product(veh[idVeh].reduced_g_prev, veh[idVeh].reduced_g_prev, numsteps, nu);
	}
	beta = nom / denom;

	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < nu; i++) {
			for (j = 0; j < numsteps; j++) {
				veh[idVeh].s_dir[i][j] = -1.0 * veh[idVeh].reduced_g[i][j] + beta * veh[idVeh].s_dir[i][j];
			}
		}
	}
}

static void polak_ribiere(struct vehicle* veh, int numsteps, int nu) {
	int i, j, idVeh;
	double beta;

	double nom = 0.0, denom = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < nu; i++) {
			for (j = 0; j < numsteps; j++) {
				veh[idVeh].grad_dif[i][j] = veh[idVeh].reduced_g[i][j] - veh[idVeh].reduced_g_prev[i][j];
			}
		}
	}

	for (idVeh = 0; idVeh < P.n; idVeh++) {
		nom += inner_product(veh[idVeh].grad_dif, veh[idVeh].reduced_g, numsteps, nu);
		denom += inner_product(veh[idVeh].reduced_g_prev, veh[idVeh].reduced_g_prev, numsteps, nu);
	}
	beta = nom / denom;

	// beta = inner_product(grad_dif, reduced_g, numsteps, NU) / inner_product(reduced_g_prev, reduced_g_prev, numsteps, NU);
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < nu; i++) {
			for (j = 0; j < numsteps; j++) {
				veh[idVeh].s_dir[i][j] = -1.0 * veh[idVeh].reduced_g[i][j] + beta * veh[idVeh].s_dir[i][j];
			}
		}
	}

}

static void rporp(struct vehicle* veh, int numsteps, int nu) {

	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		for (int i = 0; i < NU; i++) {
			for (int j = 0; j < numsteps; j++) {

				if (veh[idVeh].reduced_g_prev[i][j] * veh[idVeh].reduced_g[i][j] > 0.0) {
					veh[idVeh].delta_u[i][j] = fmin(AP.eta_plus * (veh[idVeh].delta_u[i][j]), AP.delta_u_max);
				}

				else if (veh[idVeh].reduced_g_prev[i][j] * veh[idVeh].reduced_g[i][j] < 0.0) {
					veh[idVeh].delta_u[i][j] = fmax(AP.eta_minus * (veh[idVeh].delta_u[i][j]), AP.delta_u_min);
				}
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

static void optimization(struct vehicle* veh, const char* criterion, int numsteps, double T, int iterations, double accuracy, double a, double sigma, double rho, double tau, int lin_it, int restart_it) {

	int i, j, it, idVeh, k;
	double cost, temp_cost, alpha, cost_prev;
	double beta = 0.0;
	double rprop_best = 1e4, rprop_c = 0.1;
	int rprop_theta = 0;

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

	cost = cost_prev = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++)
		cost += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

	for (idVeh = 0; idVeh < P.n; idVeh++) {
		for (i = 0; i < NU; i++) {
			for (j = 0; j < numsteps; j++) {
				veh[idVeh].reduced_g_prev[i][j] = veh[idVeh].reduced_g[i][j];
			}
		}
	}

	grad_norm[0] = grad_norm[1] = 0.0;
	for (idVeh = 0; idVeh < P.n; idVeh++)
		inner_product_1(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, NU, grad_norm);
	fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", 0, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);

	for (it = 0; it < iterations; it++) {
		if ((it == 0 || r.apply || it % restart_it == 0) && (SDIR != 5)) {
			fprintf(stderr, "restart 1 \n");
			/* Steepest Descent */
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (j = 0; j < numsteps; j++) {
						veh[idVeh].s_dir[i][j] = -1.0 * veh[idVeh].reduced_g[i][j];
					}
				}
			}
		}
		else {
			/* Calculate search direction based on the chosen algorithm */
			switch (SDIR) {
			case 1:
				/* Steepest Descent */
				for (idVeh = 0; idVeh < P.n; idVeh++) {
					for (i = 0; i < NU; i++) {
						for (j = 0; j < numsteps; j++) {
							veh[idVeh].s_dir[i][j] = -1.0 * veh[idVeh].reduced_g[i][j];
						}
					}
				}
				break;
			case 2:
				/* Fletcher - Reeves */
				fletcher_reeves(veh, numsteps, NU);
				break;
			case 3:
				/* Polak - Ribiere */
				polak_ribiere(veh, numsteps, NU);
				break;
			case 4:
				/* DFP */
				polak_ribiere(veh, numsteps, NU);
				break;
			case 5:
				/* RPROP */
				rporp(veh, numsteps, NU);
				break;
			default:
				printf("/* Invalid search direction method */\n");
				return;
			}
		}

		/* if the restart condition is satisfied break and start a new iteration */
		/* with Descent Steepest algorithm */
		if (SDIR != 5) {
			r.apply = false;
			if (it != 0 && restart(veh, numsteps, NU, it)) {
				r.apply = true;
				fprintf(stderr, "restart 2 \n");
				continue;
			}

			/* Calculate the optimum scalar step (alpha) through Line Optimization */
			alpha = line_search(veh, a, sigma, tau, rho, lin_it, numsteps, T, it, cost_prev);

			/* restart if line search limit reached */
			if (r.apply) {
				fprintf(stderr, "--restart due to line search limit reached \n");
				continue;
			}

			/* Calculate the new control */
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (j = 0; j < numsteps; j++) {
						/*if (i == 1)
							fprintf(stderr, "(%d, %d) alpha: %.15f \t sdir: %.15f \t gradient: %.15f \n", idVeh, j, alpha, veh[idVeh].s_dir[i][j], veh[idVeh].reduced_g[i][j]);*/
						veh[idVeh].temp_control[i][j] = veh[idVeh].control[i][j] + alpha * veh[idVeh].s_dir[i][j];
					}
				}
			}
		}

		for (idVeh = 0; idVeh < P.n; idVeh++)
			control_bound(veh[idVeh].temp_control, c_min, c_max, numsteps, NU, idVeh);

		/* Calculate the new states and cost */
		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].temp_control, numsteps, T, veh[idVeh].temp_states, idVeh);
		/*if (it == 215) {
			for (i = 0; i < P.numsteps; i++) {
				fprintf(stderr, "(%d) uxLB: %.4f \t uxUB: %.4f \t uyLB: %.4f \t uyUB; %.4f \n", i, veh[0].uxLB[i], veh[0].uxUB[i], veh[0].uyLB[i], veh[0].uyUB[i]);
			}
		}	*/

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
		for (idVeh = 0; idVeh < P.n; idVeh++)
			temp_cost += compute_cost(veh, states_all_veh, control_all_veh, numsteps, T, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_costates(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].temp_costates, veh[idVeh].temp_reduced_g, idVeh);

		for (idVeh = 0; idVeh < P.n; idVeh++)
			compute_reduced_gradient(veh, states_all_veh, control_all_veh, numsteps, T, veh[idVeh].temp_costates, veh[idVeh].temp_reduced_g, idVeh);

		if (termination_1(veh, criterion, accuracy, temp_cost, cost, numsteps, grad_norm, it) || it == iterations - 1) {

			// fprintf(stderr, "termination \n");
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NX; i++) {
					for (j = 0; j < numsteps + 1; j++) {
						veh[idVeh].states[i][j] = veh[idVeh].temp_states[i][j];
					}
				}
			}
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (j = 0; j < numsteps; j++) {
						veh[idVeh].control[i][j] = veh[idVeh].temp_control[i][j];
					}
				}
			}
			cost = temp_cost;

			//fprintf(stderr, "(CPU Time: %.4f)\n", ((double)(clock() - start_opt)) / CLOCKS_PER_SEC);
			fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", it + 1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);

			break;

		}
		else {
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NU; i++) {
					for (j = 0; j < numsteps; j++) {
						veh[idVeh].reduced_g_prev[i][j] = veh[idVeh].reduced_g[i][j];
						veh[idVeh].reduced_g[i][j] = veh[idVeh].temp_reduced_g[i][j];
						/* added rprop */
						veh[idVeh].control_prev[i][j] = veh[idVeh].control[i][j];
						veh[idVeh].control[i][j] = veh[idVeh].temp_control[i][j];
					}
				}
			}
			for (idVeh = 0; idVeh < P.n; idVeh++) {
				for (i = 0; i < NX; i++) {
					for (j = 0; j < numsteps + 1; j++) {
						veh[idVeh].states[i][j] = veh[idVeh].temp_states[i][j];
					}
				}
			}
			cost_prev = cost;
			cost = temp_cost;
		}

		grad_norm[0] = grad_norm[1] = 0.0;
		for (idVeh = 0; idVeh < P.n; idVeh++)
			inner_product_1(veh[idVeh].reduced_g, veh[idVeh].reduced_g, numsteps, NU, grad_norm);

		//fprintf(stderr, "(CPU Time: %.4f)\n", ((double)(clock() - start_opt)) / CLOCKS_PER_SEC);
		fprintf(stderr, "iter: %d \t cost: %.7f \t \t gradient: %.4f -- %.4f \t (CPU Time: %.4f) \t brac_iter: %d \t sec_iter: %d\n", it + 1, cost, grad_norm[0], grad_norm[1], ((double)(clock() - start_opt)) / CLOCKS_PER_SEC, brac_it, sect_it);
	}
}

static void allocations(int numsteps) {
	int i, j, idVeh;


	((veh = (struct vehicle*)calloc(P.n, sizeof(struct vehicle))));
	for (idVeh = 0; idVeh < P.n; idVeh++) {
		((veh[idVeh].control = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_prev = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].s_dir = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].temp_control = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].temp_reduced_g = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].control_prev = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].control_cur = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_cur = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_a1 = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_a2 = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].reduced_g_alpha = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].control_a1 = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].control_a2 = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].control_alpha = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].grad_dif = (double**)calloc(NU, sizeof(double*))));

		((veh[idVeh].delta_u = (double**)calloc(NU, sizeof(double*))));
		((veh[idVeh].rprop_best_control = (double**)calloc(NU, sizeof(double*))));

		for (i = 0; i < NU; i++) {
			((veh[idVeh].control[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_prev[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].s_dir[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].temp_control[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].temp_reduced_g[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].control_prev[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].control_cur[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_cur[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_a1[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_a2[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].reduced_g_alpha[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].control_a1[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].control_a2[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].control_alpha[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].grad_dif[i] = (double*)calloc(numsteps, sizeof(double))));

			((veh[idVeh].delta_u[i] = (double*)calloc(numsteps, sizeof(double))));
			((veh[idVeh].rprop_best_control[i] = (double*)calloc(numsteps, sizeof(double))));

		}

		((veh[idVeh].uxLB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uxUB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uyLB = (double*)calloc(numsteps, sizeof(double))));
		((veh[idVeh].uyUB = (double*)calloc(numsteps, sizeof(double))));

		((veh[idVeh].states = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].temp_states = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].temp_costates = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].states_prev = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].states_cur = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates_cur = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].states_a1 = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].states_a2 = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].states_alpha = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates_a1 = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates_a2 = (double**)calloc(NX, sizeof(double*))));
		((veh[idVeh].costates_alpha = (double**)calloc(NX, sizeof(double*))));
		for (i = 0; i < NX; i++) {
			((veh[idVeh].states[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].temp_states[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].temp_costates[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].states_prev[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].states_cur[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates_cur[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].states_a1[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].states_a2[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].states_alpha[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates_a1[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates_a2[i] = (double*)calloc(numsteps + 1, sizeof(double))));
			((veh[idVeh].costates_alpha[i] = (double*)calloc(numsteps + 1, sizeof(double))));
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
	double dval;
	int ival;

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

	FILE* input;
	input = fopen("input.txt", "r");
	//fopen_s(&input, "../input-0.txt", "r");
	// fopen_s(&input, "../input-1.txt", "r");

	while (fgets(buf, sizeof(buf), input))
	{
		// fprintf(stderr, "parsing: %s\n", buf);
		if (buf[0] == '\n')
			break;

		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			P.numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			P.T = dval;
		if (sscanf(buf, "\"roadwidth\":%lf", &dval) == 1)
			P.roadwidth = dval;
		if (sscanf(buf, "\"safety\":%lf", &dval) == 1)
			P.safety = dval;
		if (sscanf(buf, "\"latmar\":%lf", &dval) == 1)
			P.latmar = dval;
		if (sscanf(buf, "\"penalty(%d)\":%lf", &ival, &dval) == 2)
			penalty[ival] = dval;
		if (sscanf(buf, "\"cmux\":%lf", &dval) == 1)
			P.cmux = dval;
		if (sscanf(buf, "\"cmuy\":%lf", &dval) == 1)
			P.cmuy = dval;
		if (sscanf(buf, "\"radius\":%d", &ival) == 1)
			P.rad = ival;
		if (sscanf(buf, "\"elc\":%lf", &dval) == 1)
			P.elc = dval;
		if (sscanf(buf, "\"slope\":%lf", &dval) == 1)
			P.slp = dval;

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
		if (sscanf(buf, "\"obst_l(%d)\":%lf", &ival, &dval) == 2) {
			P.obst_l[ival] = dval;
			P.obst_l_len = MAX(P.obst_l_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_w(%d)\":%lf", &ival, &dval) == 2) {
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
		if (sscanf(buf, "\"FDA_delta_0\":%lf", &dval) == 1)
			AP.delta_u_0 = dval;
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

static void printsolSim(struct vehicle* veh, int numsteps, double T) {

	FILE* fd2 = NULL;
	if ((fd2 = fopen("../sim/data/sim.js", "w")) == NULL);

	fputs("sim = {\n", fd2);
	fputs("x: [\n", fd2);
	double max_x = 0.0;
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++) {
			fprintf(fd2, "%.4f,", veh[idVeh].states[0][k]);
			max_x = MAX(max_x, veh[idVeh].states[0][k]);
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
			fprintf(fd2, "%.8f,", veh[idVeh].control[0][k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("ay: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.8f,", veh[idVeh].control[1][k]);
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

	fputs("vy: [\n", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < P.numsteps; k++)
			fprintf(fd2, "%.4f,", veh[idVeh].states[3][k]);
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
		fprintf(fd2, "%.4f,", 5.0);
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < P.n; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	for (int idVeh = 0; idVeh < P.n_obs; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", P.n + P.n_obs);
	fprintf(fd2, "k:%d,\n", P.numsteps);
	fprintf(fd2, "roadlength:%f,\n", 5000.0);
	fprintf(fd2, "roadwidth:%f,\n", P.roadwidth);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", P.T);
	fprintf(fd2, "roadbound:%f,\n", 1.2);

	fprintf(fd2, "offset_count:%d,\n", 0);
	fprintf(fd2, "limits: [],\n");
	fprintf(fd2, "slopes: [],\n");
	fprintf(fd2, "offsets: [],\n");
	fprintf(fd2, "};\n\n");

	fclose(fd2);
}

int main() {

	read_P();
	
	P.uy_min = -1.0;
	P.uy_max =  1.0;

	P.ux_min = -2.0;
	P.ux_max =  2.0;

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
					veh[idVeh].delta_u[i][j] = AP.delta_u_0;
				}
			}
		}
	}

	double lm = 0, yLB, yUB, vyLB, vyUB;
	for (int idVeh = 0; idVeh < P.n; idVeh++) {
		lm = 2.0 / 2.0;
		fprintf(stderr, "lm: %.4f \n", lm);
		veh[idVeh].yLim[0][0] = lm;
		veh[idVeh].yLim[1][0] = P.roadwidth - lm;
		veh[idVeh].vyLim[0][0] = veh[idVeh].vyLim[1][0] = 0.0;

		for (int i = 0; i < P.numsteps; i++) {
			veh[idVeh].vyLim[0][i + 1] = vyLB = veh[idVeh].vyLim[0][0];
			veh[idVeh].vyLim[1][i + 1] = vyUB = veh[idVeh].vyLim[1][0];
			veh[idVeh].yLim[0][i + 1] = yLB = veh[idVeh].yLim[0][0];
			veh[idVeh].yLim[1][i + 1] = yUB = veh[idVeh].yLim[1][0];

			veh[idVeh].control[0][i] = (fmax(-LONK3 * P.vx0[idVeh], P.ux_min) / (-P.ux_max + fmax(-LONK3 * P.vx0[idVeh], P.ux_min)));
			veh[idVeh].control[1][i] = ((LATK1 * (-P.y0[idVeh] + yLB) + LATK2 * (vyLB - P.vy0[idVeh]))
				/ ((LATK1 * (-P.y0[idVeh] + yLB) + LATK2 * (vyLB - P.vy0[idVeh])) - (LATK1 * (-P.y0[idVeh] + yUB) + LATK2 * (vyUB - P.vy0[idVeh]))));

		}
	}

	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_states(veh, veh[idVeh].states, states_all_veh, veh[idVeh].control, P.numsteps, P.T, veh[idVeh].states, idVeh); }

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

	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_costates(veh, states_all_veh, control_all_veh, P.numsteps, P.T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh); }
	for (int idVeh = 0; idVeh < P.n; idVeh++) { compute_reduced_gradient(veh, states_all_veh, control_all_veh, P.numsteps, P.T, veh[idVeh].costates, veh[idVeh].reduced_g, idVeh); }

	start_FDA = clock();
	optimization(veh, criterion, P.numsteps, P.T, AP.iter, AP.accur, AP.a, AP.sigma, AP.rho, AP.tau, AP.lin_it, AP.restart_it);
	end_FDA = clock();
	cpu_time_used_FDA = ((double)(end_FDA - start_FDA)) / CLOCKS_PER_SEC;
	fprintf(stderr, "CPU Time: %.4f\n", cpu_time_used_FDA);

	printsolSim(veh, P.numsteps, P.T);

	return 0;
}