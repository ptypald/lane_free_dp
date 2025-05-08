#include <algorithm>
#include <float.h>
#include "dddp.h"

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

void DDDP::feasible_domain(int it) {
	vector<double> temp;

	// clear vectors
	this->UX.clear(); this->UY.clear();
	this->vxmin.clear(); this->vxmax.clear(); this->NVX.clear(); this->VX.clear();
	this->xmin.clear(); this->xmax.clear(); this->NX.clear(); this->X.clear();
	this->ymin.clear(); this->ymax.clear(); this->NY.clear(); this->Y.clear();

	/* ux and uy domain */
	this->dux = this->dux0;
	this->uxmin = uxLB; this->uxmax = uxUB;
	this->duy = this->duy0;
	this->uymin = uyLB; this->uymax = uyUB;

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
		if (fmod(init_vx[numsteps-1], this->dvx) != 0) {
			this->vxmin.push_back(floor(MAX(init_vx[k] - this->cvx * this->dvx, vxLB)));
			this->vxmax.push_back(ceil(MIN(init_vx[k] + this->cvx * this->dvx, vxUB)));
		}
		else {
			this->vxmin.push_back(MAX(init_vx[k] - this->cvx * this->dvx, vxLB));
			this->vxmax.push_back(MIN(init_vx[k] + this->cvx * this->dvx, vxUB));
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
	dx = 0.5 * dux * step;
	cx = cvx0 * cxm;
	for (int k = 0; k < numsteps + 1; k++) {
		if (fmod(init_x[numsteps-1], dx) != 0) {
			// xmin.push_back(floor(MAX(init_x[k] - (cx * dx), roadxMin)));
			// xmax.push_back(ceil(MIN(init_x[k] + (cx * dx), roadxMax)));
			xmin.push_back(floor(MAX(init_x[k] - (cx * dx), x0)));
			xmax.push_back(ceil(MIN(init_x[k] + (cx * dx), roadxMax)));
		}
		else {
			// xmin.push_back(MAX(init_x[k] - (cx * dx), roadxMin));
			// xmax.push_back(MIN(init_x[k] + (cx * dx), roadxMax));
			xmin.push_back(MAX(init_x[k] - (cx * dx), x0));
			xmax.push_back(MIN(init_x[k] + (cx * dx), roadxMax));
		}
		NX.push_back((int)ceil((xmax[k] - xmin[k]) / dx) + 1);
		
		for (int ix = 0; ix < NX[k]; ix++) {
			temp.push_back(xmin[k] + (dx * ix));
		}
		X.push_back(temp);
		temp.clear();
	}
	maxNX = *max_element(NX.begin(), NX.end());

	/* y domain */
	this->dy = 0.5 * this->duy * this->step;
	this->cy = this->cy0;
	for (int k = 0; k < this->numsteps + 1; k++) {	
		this->ymin.push_back(MAX(init_y[k] - this->cy * this->dy, roadyMin + 1.0));
		this->ymax.push_back(MIN(init_y[k] + this->cy * this->dy, roadyMax - 1.0));

		this->NY.push_back((int)ceil((this->ymax[k] - this->ymin[k]) / this->dy) + 1);

		for (int iy = 0; iy < this->NY[k]; iy++)
			temp.push_back(this->ymin[k] + this->dy * iy);
		this->Y.push_back(temp);
		temp.clear();
	}
	maxNY = *max_element(NY.begin(), NY.end());
	// plotDomain(v, it);
}

void DDDP::initial_path() {
	double uxConst;
	init_x.push_back(x0); init_y.push_back(y0);
	init_vx.push_back(vx0);
	double T = step;
	cout << "vdx: " << vdx << endl;
	uxConst = (vdx - vx0) / numsteps;
	for (int k = 0; k < numsteps; k++) {
		init_x.push_back(init_x[k] + init_vx[k]*T + 0.5*uxConst*pow(T,2));
		init_y.push_back(init_y[k]);
		init_vx.push_back(init_vx[k] + uxConst*T);
	}
	// TODO: Implement IDM or other method for initial trajectory

	// IDM - MOBIL

	fprintf(stderr, "-- Initial path created!\n");
}

void DDDP::execute(Vehicle v) {
	int it = 0;
	double x, y, vx;

	clock_t it_tic, it_toc;
	double it_cpu_time;

	v.obst_prediction();

	while(true) {
		it_tic = clock();

		feasible_domain(it);

		// TODO:
		myAlloc();

		for (int k = numsteps; k >= 0; k--) {
			// for each state variable
			for (int ix = 0; ix < ((k == 0) ? 1 : NX[k]); ix++) {
				x = ((k==0) ? x0 : X[k][ix]);
				for (int iy = 0; iy < ((k == 0) ? 1 : NY[k]); iy++) {
					y = ((k==0) ? y0 : Y[k][iy]);
					for (int ivx = 0; ivx < ((k == 0) ? 1 : NVX[k]); ivx++) {
						vx = ((k==0) ? vx0 : VX[k][ivx]);

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
						if (k == numsteps) {
							if (collision(v, here.x, here.y, here.vx, here.k) == 1) {
								J[here.k][here.ix][here.iy][here.ivx] = DBL_MAX;
							}
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

						// TODO:
						if (collision(v, here.x, here.y, here.vx, here.k) == 1)
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

								phi = 0.5*pow(ux, 2) + pow(here.vx - vdx, 2) + pow(uy - vdy, 2);
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
		// TODO:
		printSolution(v, it);

		it_toc = clock();
		it_cpu_time = (double)(it_toc - it_tic) / CLOCKS_PER_SEC;
		fprintf(stderr, "Iteration %d | Δux=Δvx: %.4f | Δx: %.5f | NX: %d | NY: %d | NVX: %d | NUX: %d | NUY: %d | cost: %.5f | CPU time: %.4f \n", 
			it + 1, dux, dx, maxNX, maxNY, maxNVX, NUX, NUY, opt_J_it, it_cpu_time);

		// if (it >= 1 || this->dux0 < 0.125) { 
		if (this->dux0 <= 0.125 && this->duy0 <= 0.125) {
			break;
		}

		this->dux0 /= 2.0; 
		this->duy0 /= 2.0;
		// re-adjust the initial trajectories
		for (int k = 0; k < numsteps + 1; k++) {
			init_x[k] = opt_x[k];
			init_y[k] = opt_y[k];
			init_vx[k] = opt_vx[k];
		}
		// clear memory and whatever else
		// TODO:
		myFree();

		it++;
	}

};
