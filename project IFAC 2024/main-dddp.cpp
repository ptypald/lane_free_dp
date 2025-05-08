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

using namespace std;

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

class Infrastructure {
	public:
		double xMin = 0.0, xMax = 100000.0, yMin = 0.0, yMax = 10.0;
};

class Data {
	public:
		// flags
		int DISPLAY = 0;
		// input
		int numsteps = 24; // 6 seconds
		double step = 1.0;
		double safety_x, safety_y;
		double x0, y0, vx0, ax0, ay0; // is this needed?
		double vdx, vdy;
		vector<double> obsx, obsy, obsvx, obsvy, obsax, obsay;
		int obst_x_len, obst_y_len, obst_vx_len, obst_vy_len;
		int obs_n;

		// output
		vector<double> x, y, vx, ux, uy; // is this needed?

		// functions
		void input();
		void exportSimData(class Vehicle veh, class Infrastructure inf, int it);
		// void exportPlotData(Vehicle veh);
};

class Vehicle {
	
	private:
		double x0, y0, vx0;
		double len, width;
		// double tgx, tgy;	// time gaps
		double vdx, vdy;	// desired speeds
		// vehicle bounds
		double UXMIN = -2.0, UXMAX = 1.0;
		double UYMIN = -1.0, UYMAX = 1.0;
		double VXMIN = 0.0, VXMAX = 40.0;

	public:
		vector<double> init_vx, init_x, init_y;
		vector<double> x, y, vx, ux, uy;
	
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
	struct Obstacles;
	vector<Obstacles> obs;

	// class functions
	void initial_path(Data data);
	void IDM(Data data);
	void obst_prediction(Data data);
	int collision(Data data, double here_x, double here_y, double here_vx, int here_k);
};

struct Vehicle::Obstacles {
	int id;
	// double x0, y0, vx0, vy0, ax0, ay0;
	vector<double> x, y, vx, vy, ax, ay;
	double len, width;
};

class Controller {

	public:
		// general parameters
		int numsteps, max_iter = 0;
		double step;
		string method;
		// domain parameters
		vector<int> NX, NVX, NY;
		int maxNX, maxNVX, maxNY;
		vector<double> xmin, xmax, vxmin, vxmax, ymin, ymax;
		vector<vector<double>> X, VX, Y;
		vector<double> UX, UY;
		double dx, dvx, dux, dy, duy;
		double uxmin, uxmax, uymin, uymax;
		int NUX, NUY;
		// corridor parameters
		double dux0 = 1.0, duy0 = 1.0, cvx0 = 5.0, cxm = 15.0, cy0 = 15.0;
		double cx, cy, cvx;
		// optimal trajectories
		vector<double> opt_x, opt_vx, opt_y, opt_ux, opt_uy;
		// 
		double opt_J_it, opt_J_prev;

		double ****J, ****Ux, ****Uy;

	// parametrized constructor
	Controller(int K, double T, string m) {
		numsteps = K;
		step = T;
		method = m;
	}

	typedef struct node {
		int k;
		int ix, ivx, iy;
		double x, vx, y;
		double _x, _vx, _y;
	} node_t;

	typedef struct {
		double ux, uy;
	} ctrl_t;

	void feasible_domain(Vehicle v, Infrastructure r, int it);
	void dddp(Vehicle v, Infrastructure r, Data s);
	void allocations();
	void myFree();
	int discretizexv(node_t *state, int k);
	void printSolution(Vehicle v, Data data, Infrastructure inf, int it);
	void plotDomain(Vehicle v, int it);
	void plotSolution(Vehicle v, int it);

	node_t *getnext(node_t here, double ux, double uy){
		static node_t next;
		int err;

		next.k = here.k + 1;
		next.x = next._x = here.x + here.vx*step + 0.5*ux*pow(step, 2);
		next.y = next._y = here.y + uy*step;
		next.vx = next._vx = here.vx + ux*step;
		
		if ((err = discretizexv(&next, next.k))) {
			fprintf(stderr, " -- INVALID STATE: form (%.3f | %.3f | %.3f) to (%.3f | %.3f | %.3f) at time step %d to %d \n", 
				here.x, here.y, here.vx, next.x, next.y, next.vx, here.k, next.k);
			return NULL;
		}
		return &next;
	};
};

void Data::input() {

	char buf[1024];
	// const char* line;
	double dval;
	int ival;

	FILE* input;
	input = fopen("input.txt", "r");

	while (fgets(buf, sizeof(buf), input)) {
		// fprintf(stderr, "parsing: %s\n", buf);
		if (buf[0] == '\n')
			break;

		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			this->numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			this->step = dval;
		if (sscanf(buf, "\"safety_x\":%lf", &dval) == 1)
			safety_x = dval;
		if (sscanf(buf, "\"safety_y\":%lf", &dval) == 1)
			safety_y = dval;

		/* read initial states */
		if (sscanf(buf, "\"vdx\":%lf", &dval) == 1)
			vdx = dval;
		if (sscanf(buf, "\"vdy\":%lf", &dval) == 1)
			vdy = dval;
		if (sscanf(buf, "\"x(0)\":%lf", &dval) == 1)
			x0 = dval;
		if (sscanf(buf, "\"y(0)\":%lf", &dval) == 1)
			y0 = dval;
		if (sscanf(buf, "\"vx(0)\":%lf", &dval) == 1)
			vx0 = dval;
		if (sscanf(buf, "\"ax(0)\":%lf", &dval) == 1)
			ax0 = dval;
		if (sscanf(buf, "\"ay(0)\":%lf", &dval) == 1)
			ay0 = dval;
		
		/* obstacles */
		if (sscanf(buf, "\"obst_x(%d,0)\":%lf", &ival, &dval) == 2)
			obsx.push_back(dval);
		if (sscanf(buf, "\"obst_y(%d,0)\":%lf", &ival, &dval) == 2)
			obsy.push_back(dval);
		if (sscanf(buf, "\"obst_vx(%d,0)\":%lf", &ival, &dval) == 2)
			obsvx.push_back(dval);
		if (sscanf(buf, "\"obst_vy(%d,0)\":%lf", &ival, &dval) == 2)
			obsvy.push_back(dval);

		memset(buf, 0, sizeof(buf));
	}

	if (!(obsx.size() == obsy.size() && obsx.size() == obsvx.size() && obsx.size() == obsvy.size())) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	obs_n = obsx.size();

}

void Vehicle::IDM(Data data) {



}

void Vehicle::initial_path(Data data) {
	
	double uxConst;
	init_x.push_back(x0); init_y.push_back(y0);
	init_vx.push_back(vx0);
	double T = data.step;
	
	uxConst = (data.vdx - data.vx0) / data.numsteps;
	for (int k = 0; k < data.numsteps; k++) {
		init_x.push_back(init_x[k] + init_vx[k]*T + 0.5*uxConst*pow(T,2));
		init_y.push_back(init_y[k]);
		init_vx.push_back(init_vx[k] + uxConst*T);
	}
	// TODO: Implement IDM or other method for initial trajectory

	// IDM - MOBIL

	fprintf(stderr, "-- Initial path created!\n");
}

void Vehicle::obst_prediction(Data data) {

	for (int i = 0; i < data.obs_n; i++) {
		Obstacles temp_obst;
		temp_obst.x.push_back(data.obsx[i]); temp_obst.y.push_back(data.obsy[i]);
		temp_obst.vx.push_back(data.obsvx[i]); temp_obst.vy.push_back(data.obsvy[i]);
		for (int k = 0; k < data.numsteps; k++) {
			temp_obst.x.push_back(temp_obst.x[k] + temp_obst.vx[k]*data.step);
			temp_obst.y.push_back(temp_obst.y[k] + temp_obst.vy[k]*data.step);
			temp_obst.vx.push_back(temp_obst.vx[k]);
			temp_obst.vy.push_back(temp_obst.vy[k]);
		}
		obs.push_back(temp_obst);
	}

	fprintf(stderr, "-- Obstacles path created!\n");
}

int Vehicle::collision(Data data, double here_x, double here_y, double here_vx, int here_k) {
	for (int i = 0; i < data.obs_n; i++) {
		if ( (here_x < obs[i].x[here_k] + here_vx*data.safety_x && here_x > obs[i].x[here_k] - data.obsvx[i]*data.safety_x) && (here_y < obs[i].y[here_k] + 2.0 && here_y > obs[i].y[here_k] - 2.0) ) {
			return 1;
		}
	}
 
	return 0;
}

void Controller::feasible_domain(Vehicle v, Infrastructure r, int it) {
	vector<double> temp;

	// clear vectors
	this->UX.clear(); this->UY.clear();
	this->vxmin.clear(); this->vxmax.clear(); this->NVX.clear(); this->VX.clear();
	this->xmin.clear(); this->xmax.clear(); this->NX.clear(); this->X.clear();
	this->ymin.clear(); this->ymax.clear(); this->NY.clear(); this->Y.clear();

	/* ux and uy domain */
	this->dux = this->dux0;
	this->uxmin = v.getUXmin(); this->uxmax = v.getUXmax();
	this->duy = this->duy0;
	this->uymin = v.getUYmin(); this->uymax = v.getUYmax(); 

	this->NUX = (int)ceil((this->uxmax - this->uxmin) / this->dux) + 1;
	this->NUY = (int)ceil((this->uymax - this->uymin) / this->duy) + 1;
	for (int iux = 0; iux < this->NUX; iux++) {
		this->UX.push_back(this->uxmin + this->dux * iux);
	}
	for (int iuy = 0; iuy < this->NUY; iuy++) {
		this->UY.push_back(this->uymin + this->duy * iuy);
	}

	/* vx domain */
	this->dvx = this->dux * this->step; this->cvx = this->cvx0;
	for (int k = 0; k < this->numsteps + 1; k++) {
		if (it == 0 && fmod(v.init_vx[numsteps-1], this->dvx) != 0) {
			this->vxmin.push_back(floor(MAX(v.init_vx[k] - this->cvx * this->dvx, v.getVXmin())));
			this->vxmax.push_back(ceil(MIN(v.init_vx[k] + this->cvx * this->dvx, v.getVXmax())));
		}
		else {
			this->vxmin.push_back(MAX(v.init_vx[k] - this->cvx * this->dvx, v.getVXmin()));
			this->vxmax.push_back(MIN(v.init_vx[k] + this->cvx * this->dvx, v.getVXmax()));
		}
		this->NVX.push_back((int)ceil((this->vxmax[k] - this->vxmin[k]) / this->dvx) + 1);
		
		for (int ivx = 0; ivx < this->NVX[k]; ivx++) {
			temp.push_back(this->vxmin[k] + this->dvx * ivx);
		}
		this->VX.push_back(temp);
		temp.clear();
	}
	maxNVX = *max_element(NVX.begin(), NVX.end());

	/* x domain */
	this->dx = 0.5 * this->dux * this->step;
	this->cx = this->cvx0 * this->cxm;
	for (int k = 0; k < this->numsteps + 1; k++) {
		if (it == 0 && fmod(v.init_x[numsteps-1], this->dx) != 0) {
			this->xmin.push_back(floor(MAX(v.init_x[k] - (this->cx * this->dx), r.xMin)));
			this->xmax.push_back(ceil(MIN(v.init_x[k] + (this->cx * this->dx), r.xMax)));
		}
		else {
			this->xmin.push_back(MAX(v.init_x[k] - (this->cx * this->dx), r.xMin));
			this->xmax.push_back(MIN(v.init_x[k] + (this->cx * this->dx), r.xMax));
		}
		this->NX.push_back((int)ceil((this->xmax[k] - this->xmin[k]) / this->dx) + 1);
		
		for (int ix = 0; ix < this->NX[k]; ix++) {
			temp.push_back(this->xmin[k] + (this->dx * ix));
		}
		this->X.push_back(temp);
		temp.clear();
	}
	maxNX = *max_element(NX.begin(), NX.end());

	/* y domain */
	this->dy = 0.5 * this->duy * this->step;
	this->cy = this->cy0;
	for (int k = 0; k < this->numsteps + 1; k++) {
		if (it == 0 && fmod(v.init_y[numsteps-1], this->dy) != 0) {
			this->ymin.push_back(floor(MAX(v.init_y[k] - (this->cy * this->dy), r.yMin)));
			this->ymax.push_back(ceil(MIN(v.init_y[k] + (this->cy * this->dy), r.yMax)));
		}
		else {
			this->ymin.push_back(MAX(v.init_y[k] - this->cy * this->dy, r.yMin + 1.0));
			this->ymax.push_back(MIN(v.init_y[k] + this->cy * this->dy, r.yMax - 1.0));
		}
		this->NY.push_back((int)ceil((this->ymax[k] - this->ymin[k]) / this->dy) + 1);

		for (int iy = 0; iy < this->NY[k]; iy++)
			temp.push_back(this->ymin[k] + this->dy * iy);
		this->Y.push_back(temp);
		temp.clear();
	}
	maxNY = *max_element(NY.begin(), NY.end());


	plotDomain(v, it);
}

int Controller::discretizexv(node_t *state, int k) {
	int ix, iy, ivx;

	ix = (int)round((state->x - this->X[k][0]) / this->dx);
	iy = (int)round((state->y - this->Y[k][0]) / this->dy);
	ivx = (int)round((state->vx - this->VX[k][0]) / this->dvx);
	
	if (!(ix >= 0 && ix < this->NX[k])) 
		return 1;
	if (!(iy >= 0 && iy < this->NY[k])) 
		return 2;
	if (!(ivx >= 0 && ivx < this->NVX[k]))
		return 3;

	state->x = this->X[k][(state->ix = ix)];
	state->y = this->Y[k][(state->iy = iy)];
	state->vx = this->VX[k][(state->ivx = ivx)];

	return 0;
}

void Controller::allocations() {

	this->J = (double ****)calloc(sizeof(double***), this->numsteps + 1);
	this->Ux = (double ****)calloc(sizeof(double***), this->numsteps + 1);
	this->Uy = (double ****)calloc(sizeof(double***), this->numsteps + 1);
	for (int k = 0; k < this->numsteps + 1; k++) {
		this->J[k] = (double ***)calloc(sizeof(double**), maxNX);
		this->Ux[k] = (double ***)calloc(sizeof(double**), maxNX);
		this->Uy[k] = (double ***)calloc(sizeof(double**), maxNX);
		for (int ix = 0; ix < this->maxNX; ix++) {
			this->J[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			this->Ux[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			this->Uy[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			for (int iy = 0; iy < this->maxNY; iy++) {
				this->J[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
				this->Ux[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
				this->Uy[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
			}
		}
	}
	return;
}

void Controller::myFree() {
	
	for (int k = 0; k < numsteps; k++) {
		for (int ix = 0; ix < maxNX; ix++) {
			for (int iy = 0; iy < maxNY; iy++) {
				free(J[k][ix][iy]);
				free(Ux[k][ix][iy]);
				free(Uy[k][ix][iy]);
			}
		}
	}

	for (int k = 0; k < numsteps; k++) {
		for (int ix = 0; ix < maxNX; ix++) {
			free(J[k][ix]);
			free(Ux[k][ix]);
			free(Uy[k][ix]);
		}
	}

	for (int k = 0; k < numsteps; k++) {
		free(J[k]);
		free(Ux[k]);
		free(Uy[k]);
	}

	free(J);
	free(Ux);
	free(Uy);


}

void Controller::printSolution(Vehicle v, Data data, Infrastructure inf, int it){

	node_t *nextptr;
	node_t *state;
	double ux, uy, real_cost = 0.0;

	data.ux.clear(); data.uy.clear();
	data.x.clear(); data.y.clear(); data.vx.clear();
	opt_ux.clear(); opt_uy.clear();
	opt_x.clear(); opt_y.clear(); opt_vx.clear();

	node_t initial = { 0 };
	initial.k = 0;
	initial.x = v.getX0();
	initial.y = v.getY0();
	initial.vx = v.getVX0();
	initial.ix = initial.iy = initial.ivx = 0;

	state = (node_t*)calloc(sizeof(node_t), this->numsteps + 1);
	state[0] = initial;

	if (data.DISPLAY) {
		fprintf(stderr, "k \t x \t\t y \t\t vx \t\t ux \t\t uy \t\t J \t\t real cost \n");
		fprintf(stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	}
	for (int k = 0; k <= this->numsteps; k++){
		ux = ACCESS(Ux, state[k]); data.ux.push_back(ux); opt_ux.push_back(ux);
		uy = ACCESS(Uy, state[k]); data.uy.push_back(uy); opt_uy.push_back(uy);
		data.x.push_back(state[k].x); data.y.push_back(state[k].y); data.vx.push_back(state[k].vx);
		opt_x.push_back(state[k].x); opt_y.push_back(state[k].y); opt_vx.push_back(state[k].vx);

		// real cost
		real_cost += 0.5*pow(ux, 2) + 0.5*pow(state[k].vx - v.getDersiredSpeedX(), 2) + 0.5*pow(uy - v.getDersiredSpeedY(), 2);
		
		if (data.DISPLAY) {
			fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
				k,
				state[k].x,
				state[k].y,
				state[k].vx,
				ux,
				uy,
				ACCESS(J, state[k]),
				real_cost
			);
		}

		if (k < this->numsteps) {
			if ((nextptr = getnext(state[k], ux, uy)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}
	}

	opt_J_it = ACCESS(J, state[0]);
	free(state);
	data.exportSimData(v, inf, it);

}

void Controller::dddp(Vehicle v, Infrastructure r, Data data) {
	int it = 0;
	double x, y, vx;
	double x0 = v.getX0(), y0 = v.getY0(), vx0 = v.getVX0();
	v.setDersiredSpeedX(data.vdx); v.setDersiredSpeedY(data.vdy);

	clock_t it_tic, it_toc;
	double it_cpu_time;

	while(true) {
		it_tic = clock();
		// TODO: check this
		this->feasible_domain(v, r, it);
		this->allocations();
		// fprintf(stderr, "(%d) dux: %.2f -- duy: %.2f\n", it+1, this->dux, this->duy);

		for (int k = this->numsteps; k >= 0; k--) {
			// for each state variable
			for (int ix = 0; ix < ((k == 0) ? 1 : this->NX[k]); ix++) {
				x = ((k==0) ? x0 : this->X[k][ix]);
				for (int iy = 0; iy < ((k == 0) ? 1 : this->NY[k]); iy++) {
					y = ((k==0) ? y0 : this->Y[k][iy]);
					for (int ivx = 0; ivx < ((k == 0) ? 1 : this->NVX[k]); ivx++) {
						vx = ((k==0) ? vx0 : this->VX[k][ivx]);
						
						// initialize current node
						node_t here = { 0 };
						here.k = k;
						here.x = x; here.y = y;
						here.vx = vx;
						here.ix = ((k == 0) ? 0 : ix); here.iy = ((k == 0) ? 0 : iy); 
						here.ivx = ((k == 0) ? 0 : ivx);
						
						if (here.ix > this->NX[k] || here.iy > this->NY[k] || here.ivx > this->NVX[k]) {
							fprintf(stderr, "invalid state on the grid\n"); 
							exit(1);
						}

						// final time step costs are zero
						if (k == this->numsteps) {
							if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
								J[here.k][here.ix][here.iy][here.ivx] = DBL_MAX;
							else
								J[here.k][here.ix][here.iy][here.ivx] = 0.0;			

							Ux[here.k][here.ix][here.iy][here.ivx] = 0.0;
							Uy[here.k][here.ix][here.iy][here.ivx] = 0.0;
							
							continue;
						}

						node_t next = { 0 };
						double cost, phi;
						double *Jhere, *Jnext, *Uxhere, *Uyhere;
						
						Jhere = &J[here.k][here.ix][here.iy][here.ivx];
						Uxhere = &Ux[here.k][here.ix][here.iy][here.ivx];
						Uyhere = &Uy[here.k][here.ix][here.iy][here.ivx];
						
						*Jhere = DBL_MAX;
						next.k = here.k + 1;

						// if (crash)
						if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
							continue;
						
						// for each control variable
						for (auto ux : UX) {
							for (auto uy : UY) {
								next.x = here.x + here.vx*this->step + 0.5*ux*pow(this->step, 2);
								next.y = here.y + uy*this->step;
								next.vx = here.vx + ux*this->step;

								// if ( feasible discrete point )
								if (discretizexv(&next, next.k) != 0)
									continue;

								// TODO: Interpolation
								
								// compute cost
								Jnext = &J[next.k][next.ix][next.iy][next.ivx];

								phi = 0.5*pow(ux, 2) + pow(here.vx - v.getDersiredSpeedX(), 2) + pow(uy - v.getDersiredSpeedY(), 2);
								cost =  phi + *Jnext;

								if (cost < *Jhere) {
									*Jhere = cost;
									*Uxhere = ux;
									*Uyhere = uy;
								}
							}					
						}
					}
				}
			}
		}
		this->printSolution(v, data, r, it);
		this->plotSolution(v, it);
		it_toc = clock();
		it_cpu_time = (double)(it_toc - it_tic) / CLOCKS_PER_SEC;
		fprintf(stderr, "Iteration %d | Δux=Δvx: %.4f | Δx: %.5f | NX: %d | NY: %d | NVX: %d | NUX: %d | NUY: %d | cost: %.5f | CPU time: %.4f \n", 
			it + 1, dux, dx, maxNX, maxNY, maxNVX, NUX, NUY, opt_J_it, it_cpu_time);

		// if (it >= 1 || this->dux0 < 0.125) { 
		// if (this->dux0 <= 0.125 && this->duy0 <= 0.125) {
		cout << "costs: " << opt_J_it << " | " << opt_J_prev << endl;
		if (fabs(opt_J_it - opt_J_prev) < 0.1) {
			break;
		}

		opt_J_prev = opt_J_it;
		this->dux0 /= 2.0; 
		this->duy0 /= 2.0;
		// re-adjust the initial trajectories
		for (int k = 0; k < numsteps + 1; k++) {
			v.init_x[k] = opt_x[k];
			v.init_y[k] = opt_y[k];
			v.init_vx[k] = opt_vx[k];
		}
		// clear memory and whatever else
		myFree();

		it++;
	}
}

void Data::exportSimData(class Vehicle veh, class Infrastructure inf, int it) {
	
	// char fileBuffer[32];
    // snprintf(fileBuffer, sizeof(char) * 32, "../../sim/data/sim_%i.js", it);

	FILE* fd2 = NULL;
	if ((fd2 = fopen("../../sim/data/sim.js", "w")) == NULL) {
		printf("\nCouldn't read data file \n");
	}

	double max_pos = 0.0;
	fprintf(fd2, "sim = {\n");
	fprintf(fd2, "x: [\n");
	
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++) {
			fprintf(fd2, "%.4f,", x[k]);
			if (x[k] > max_pos)
				max_pos = x[k];
		}
		fputs("],\n", fd2);
	}
	
	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++) {
			fprintf(fd2, "%.4f,", veh.obs[idVeh].x[k]);
			if (veh.obs[idVeh].x[k] > max_pos)
				max_pos = veh.obs[idVeh].x[k];
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("ax: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", this->ux[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("ay: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", this->uy[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("y: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", this->y[k]);
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", veh.obs[idVeh].y[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vx: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", this->vx[k]);
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < this->numsteps; k++)
			fprintf(fd2, "%.4f,", veh.obs[idVeh].vx[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vdx: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%.4f,", veh.getDersiredSpeedX());
	fputs("],\n", fd2);

	fputs("id: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%d,", idVeh + 1);
	fputs("],\n", fd2);

	fputs("Cx: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	for (int idVeh = 0; idVeh < this->obs_n; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	for (int idVeh = 0; idVeh < this->obs_n; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", 1 + this->obs_n);
	fprintf(fd2, "k:%d,\n", this->numsteps);
	fprintf(fd2, "roadlength:%f,\n", max_pos);
	fprintf(fd2, "roadwidth:%f,\n", inf.yMax - inf.xMin);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", this->step);
	fprintf(fd2, "safety1:%f,\n", this->safety_x);
	fprintf(fd2, "safety2:%f,\n", this->safety_y);
	fprintf(fd2, "roadbound:%f,\n", 1.2);

	fprintf(fd2, "};\n\n");

	fclose(fd2);
}

void Controller::plotDomain(Vehicle v, int it) {
	char fileBuffer[32];
    snprintf(fileBuffer, sizeof(char) * 32, "outputs/dddp/domain/domain_%i.py", it);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");

	/* python file */
	fputs("domain = {\n", fout);
	
	// initial x
	fputs("'x0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", v.init_x[k]);
	fputs("],\n", fout);

	// initial y
	fputs("'y0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", v.init_y[k]);
	fputs("],\n", fout);

	// initial vx
	fputs("'vx0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", v.init_vx[k]);
	fputs("],\n", fout);
	
	// vx
	fputs("'vxmin': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", vxmin[k]);
	fputs("],\n", fout);

	fputs("'vxmax': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", vxmax[k]);
	fputs("],\n", fout);

	fputs("'vx': [\n", fout);
	for (int k = 0; k < numsteps; k++) {
		fputs("\t[", fout);
		for (int ivx = 0; ivx < NVX[k]; ivx++)
			fprintf(fout, "%.4f,", VX[k][ivx]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	// x
	fputs("'xmin': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", xmin[k]);
	fputs("],\n", fout);

	fputs("'xmax': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", xmax[k]);
	fputs("],\n", fout);

	fputs("'x': [\n", fout);
	for (int k = 0; k < numsteps; k++) {
		fputs("\t[", fout);
		for (int ivx = 0; ivx < NX[k]; ivx++)
			fprintf(fout, "%.4f,", X[k][ivx]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	// y
	fputs("'ymin': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", ymin[k]);
	fputs("],\n", fout);

	fputs("'ymax': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", ymax[k]);
	fputs("],\n", fout);

	fputs("'y': [\n", fout);
	for (int k = 0; k < numsteps; k++) {
		fputs("\t[", fout);
		for (int ivx = 0; ivx < NY[k]; ivx++)
			fprintf(fout, "%.4f,", Y[k][ivx]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fprintf(fout, "'k':%d,\n", numsteps);
	fprintf(fout, "'Step':%f,\n", step);
	fprintf(fout, "'iter':%d,\n", it + 1);
	fprintf(fout, "'dux':%f,\n", dux);
	fprintf(fout, "'duy':%f,\n", duy);
	fprintf(fout, "'dvx':%f,\n", dvx);
	fprintf(fout, "'dx':%f,\n", dx);
	fprintf(fout, "'dy':%f,\n", dy);
	fprintf(fout, "}\n\n");


	fclose(fout);
}

void Controller::plotSolution(Vehicle v, int it) {
	char fileBuffer[32];
    snprintf(fileBuffer, sizeof(char) * 32, "outputs/dddp/solution/solution_%i.py", it);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");

	/* python file */
	fputs("solution = {\n", fout);
	
	// vx
	fputs("'vx': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_vx[k]);
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

	// uy
	fputs("'uy': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_uy[k]);
	fputs("],\n", fout);

	fprintf(fout, "'k':%d,\n", numsteps);
	fprintf(fout, "'Step':%f,\n", step);
	fprintf(fout, "'iter':%d,\n", it + 1);
	fprintf(fout, "}\n\n");


	fclose(fout);
}

int main() {

	clock_t tic, toc;
	double cpu_time;
	
	tic = clock();
	/* start procedure */

	Data data;
	data.input();

	Infrastructure road;

	Vehicle veh(data.x0, data.y0, data.vx0);
	veh.obst_prediction(data);
	veh.initial_path(data);

	Controller ctr(data.numsteps, data.step, "dddp");
	ctr.dddp(veh, road, data);

	/* end procedure */
	toc = clock();
	cpu_time = (double)(toc - tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- Total CPU Time used: %.6f --\n", cpu_time);

	return 0;
}
