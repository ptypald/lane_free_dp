#include "controllers.h"

// void Controller::feasible_domain(Vehicle v, Infrastructure r, int it) {
// 	vector<double> temp;

// 	// clear vectors
// 	this->UX.clear(); this->UY.clear();
// 	this->vxmin.clear(); this->vxmax.clear(); this->NVX.clear(); this->VX.clear();
// 	this->xmin.clear(); this->xmax.clear(); this->NX.clear(); this->X.clear();
// 	this->ymin.clear(); this->ymax.clear(); this->NY.clear(); this->Y.clear();

// 	/* ux and uy domain */
// 	this->dux = this->dux0;
// 	this->uxmin = v.getUXmin(); this->uxmax = v.getUXmax();
// 	this->duy = this->duy0;
// 	this->uymin = v.getUYmin(); this->uymax = v.getUYmax(); 

// 	this->NUX = (int)ceil((this->uxmax - this->uxmin) / this->dux) + 1;
// 	this->NUY = (int)ceil((this->uymax - this->uymin) / this->duy) + 1;
// 	for (int iux = 0; iux < this->NUX; iux++) {
// 		this->UX.push_back(this->uxmin + this->dux * iux);
// 	}
// 	for (int iuy = 0; iuy < this->NUY; iuy++) {
// 		this->UY.push_back(this->uymin + this->duy * iuy);
// 	}

// 	/* vx domain */
// 	// this->dvx = this->dux * this->step; 
// 	// this->cvx = this->cvx0 * this->dvx;
// 	this->dvx = this->dux * this->step; 
// 	this->cvx = this->cvx0;
// 	for (int k = 0; k < this->numsteps + 1; k++) {
// 		if (it == 0 && fmod(v.init_vx[numsteps-1], this->dvx) != 0) {
// 			this->vxmin.push_back(floor(MAX(v.init_vx[k] - this->cvx, v.getVXmin())));
// 			this->vxmax.push_back(ceil(MIN(v.init_vx[k] + this->cvx, v.getVXmax())));
// 		}
// 		else {
// 			this->vxmin.push_back(MAX(v.init_vx[k] - this->cvx, v.getVXmin()));
// 			this->vxmax.push_back(MIN(v.init_vx[k] + this->cvx, v.getVXmax()));
// 		}
// 		this->NVX.push_back((int)ceil((this->vxmax[k] - this->vxmin[k]) / this->dvx) + 1);
		
// 		for (int ivx = 0; ivx < this->NVX[k]; ivx++) {
// 			temp.push_back(this->vxmin[k] + this->dvx * ivx);
// 		}
// 		this->VX.push_back(temp);
// 		temp.clear();
// 	}
// 	maxNVX = *max_element(NVX.begin(), NVX.end());

// 	/* x domain */
// 	// this->dx = 0.5 * this->dux * this->step;
// 	// this->cx = this->cvx0 * this->cxm * this->dx;
// 	this->dx = 0.5 * this->dux * this->step;
// 	this->cx = this->cx0;
// 	for (int k = 0; k < this->numsteps + 1; k++) {
// 		if (it == 0 && fmod(v.init_x[numsteps-1], this->dx) != 0) {
// 			this->xmin.push_back(floor(MAX(v.init_x[k] - this->cx, r.xMin)));
// 			this->xmax.push_back(ceil(MIN(v.init_x[k] + this->cx, r.xMax)));
// 		}
// 		else {
// 			this->xmin.push_back(MAX(v.init_x[k] - this->cx, r.xMin));
// 			this->xmax.push_back(MIN(v.init_x[k] + this->cx, r.xMax));
// 		}
// 		this->NX.push_back((int)ceil((this->xmax[k] - this->xmin[k]) / this->dx) + 1);
		
// 		for (int ix = 0; ix < this->NX[k]; ix++) {
// 			temp.push_back(this->xmin[k] + (this->dx * ix));
// 		}
// 		this->X.push_back(temp);
// 		temp.clear();
// 	}
// 	maxNX = *max_element(NX.begin(), NX.end());

// 	/* y domain */
// 	// this->dy = 0.5 * this->duy * this->step;
// 	// this->cy = this->cy0 * this->dy;
// 	this->dy = this->duy * this->step;
// 	this->cy = this->cy0;
// 	for (int k = 0; k < this->numsteps + 1; k++) {
// 		if (it == 0 && fmod(v.init_y[numsteps-1], this->dy) != 0) {
// 			this->ymin.push_back(floor(MAX(v.init_y[k] - this->cy, r.yMin + 1.0))); // +1 and -1 for the dimensions of the vehicle (TODO: add the real dimensions)
// 			this->ymax.push_back(ceil(MIN(v.init_y[k] + this->cy, r.yMax - 1.0)));
// 		}
// 		else {
// 			this->ymin.push_back(MAX(v.init_y[k] - this->cy, r.yMin + 1.0));
// 			this->ymax.push_back(MIN(v.init_y[k] + this->cy, r.yMax - 1.0));
// 		}
// 		this->NY.push_back((int)ceil((this->ymax[k] - this->ymin[k]) / this->dy) + 1);

// 		for (int iy = 0; iy < this->NY[k]; iy++)
// 			temp.push_back(this->ymin[k] + this->dy * iy);
// 		this->Y.push_back(temp);
// 		temp.clear();
// 	}
// 	maxNY = *max_element(NY.begin(), NY.end());


// 	plotDomain(v, it);
// }

// int Controller::discretizexv(node_t *state, int k) {
// 	int ix, iy, ivx;

// 	ix = (int)round((state->x - this->X[k][0]) / this->dx);
// 	iy = (int)round((state->y - this->Y[k][0]) / this->dy);
// 	ivx = (int)round((state->vx - this->VX[k][0]) / this->dvx);
	
// 	if (!(ix >= 0 && ix < this->NX[k])) 
// 		return 1;
// 	if (!(iy >= 0 && iy < this->NY[k])) 
// 		return 2;
// 	if (!(ivx >= 0 && ivx < this->NVX[k]))
// 		return 3;

// 	state->x = this->X[k][(state->ix = ix)];
// 	state->y = this->Y[k][(state->iy = iy)];
// 	state->vx = this->VX[k][(state->ivx = ivx)];

// 	return 0;
// }

// void Controller::allocations() {

// 	this->J = (double ****)calloc(sizeof(double***), this->numsteps + 1);
// 	this->Ux = (double ****)calloc(sizeof(double***), this->numsteps + 1);
// 	this->Uy = (double ****)calloc(sizeof(double***), this->numsteps + 1);
// 	for (int k = 0; k < this->numsteps + 1; k++) {
// 		this->J[k] = (double ***)calloc(sizeof(double**), maxNX);
// 		this->Ux[k] = (double ***)calloc(sizeof(double**), maxNX);
// 		this->Uy[k] = (double ***)calloc(sizeof(double**), maxNX);
// 		for (int ix = 0; ix < this->maxNX; ix++) {
// 			this->J[k][ix] = (double **)calloc(sizeof(double*), maxNY);
// 			this->Ux[k][ix] = (double **)calloc(sizeof(double*), maxNY);
// 			this->Uy[k][ix] = (double **)calloc(sizeof(double*), maxNY);
// 			for (int iy = 0; iy < this->maxNY; iy++) {
// 				this->J[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
// 				this->Ux[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
// 				this->Uy[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
// 			}
// 		}
// 	}
// 	return;
// }

// void Controller::myFree() {
	
// 	for (int k = 0; k < numsteps; k++) {
// 		for (int ix = 0; ix < maxNX; ix++) {
// 			for (int iy = 0; iy < maxNY; iy++) {
// 				free(J[k][ix][iy]);
// 				free(Ux[k][ix][iy]);
// 				free(Uy[k][ix][iy]);
// 			}
// 		}
// 	}

// 	for (int k = 0; k < numsteps; k++) {
// 		for (int ix = 0; ix < maxNX; ix++) {
// 			free(J[k][ix]);
// 			free(Ux[k][ix]);
// 			free(Uy[k][ix]);
// 		}
// 	}

// 	for (int k = 0; k < numsteps; k++) {
// 		free(J[k]);
// 		free(Ux[k]);
// 		free(Uy[k]);
// 	}

// 	free(J);
// 	free(Ux);
// 	free(Uy);


// }

// void Controller::printSolution(Vehicle v, Data data, Infrastructure inf, int it){

// 	node_t *nextptr;
// 	node_t *state;
// 	double ux, uy, real_cost = 0.0;

// 	data.ux.clear(); data.uy.clear();
// 	data.x.clear(); data.y.clear(); data.vx.clear();
// 	opt_ux.clear(); opt_uy.clear();
// 	opt_x.clear(); opt_y.clear(); opt_vx.clear();

// 	node_t initial = { 0 };
// 	initial.k = 0;
// 	initial.x = v.getX0();
// 	initial.y = v.getY0();
// 	initial.vx = v.getVX0();
// 	initial.ix = initial.iy = initial.ivx = 0;

// 	state = (node_t*)calloc(sizeof(node_t), this->numsteps + 1);
// 	state[0] = initial;

// 	if (data.DISPLAY) {
// 		fprintf(stderr, "k \t x \t\t y \t\t vx \t\t ux \t\t uy \t\t J \t\t real cost \n");
// 		fprintf(stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
// 	}
// 	for (int k = 0; k <= this->numsteps; k++){
// 		ux = ACCESS(Ux, state[k]); data.ux.push_back(ux); opt_ux.push_back(ux);
// 		uy = ACCESS(Uy, state[k]); data.uy.push_back(uy); opt_uy.push_back(uy);
// 		data.x.push_back(state[k].x); data.y.push_back(state[k].y); data.vx.push_back(state[k].vx);
// 		opt_x.push_back(state[k].x); opt_y.push_back(state[k].y); opt_vx.push_back(state[k].vx);

// 		// real cost
// 		real_cost += 0.5*pow(ux, 2) + 0.5*pow(state[k].vx - v.getDersiredSpeedX(), 2) + 0.5*pow(uy - v.getDersiredSpeedY(), 2);
		
// 		if (data.DISPLAY) {
// 			fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
// 				k,
// 				state[k].x,
// 				state[k].y,
// 				state[k].vx,
// 				ux,
// 				uy,
// 				ACCESS(J, state[k]),
// 				real_cost
// 			);
// 		}

// 		if (k < this->numsteps) {
// 			if ((nextptr = getnext(state[k], ux, uy)) == NULL) {
// 				fprintf(stderr, "invalid \n");
// 				exit(1);
// 			}
// 			state[k + 1] = *nextptr;
// 		}
// 	}

// 	opt_J_it = ACCESS(J, state[0]);
// 	free(state);
// 	// data.exportSimData(v, inf, it, numsteps);

// }

// void Controller::dddp(Vehicle v, Infrastructure r, Data data) {
// 	int it = 0;
// 	double x, y, vx;
// 	double x0 = v.getX0(), y0 = v.getY0(), vx0 = v.getVX0();
// 	v.setDersiredSpeedX(data.vdx); v.setDersiredSpeedY(data.vdy);

// 	clock_t it_tic, it_toc;
// 	double it_cpu_time;

// 	while(true) {
// 		it_tic = clock();
// 		// TODO: check this
// 		this->feasible_domain(v, r, it);
// 		this->allocations();
// 		// fprintf(stderr, "(%d) dux: %.2f -- duy: %.2f\n", it+1, this->dux, this->duy);

// 		for (int k = this->numsteps; k >= 0; k--) {
// 			// for each state variable
// 			for (int ix = 0; ix < ((k == 0) ? 1 : this->NX[k]); ix++) {
// 				x = ((k==0) ? x0 : this->X[k][ix]);
// 				for (int iy = 0; iy < ((k == 0) ? 1 : this->NY[k]); iy++) {
// 					y = ((k==0) ? y0 : this->Y[k][iy]);
// 					for (int ivx = 0; ivx < ((k == 0) ? 1 : this->NVX[k]); ivx++) {
// 						vx = ((k==0) ? vx0 : this->VX[k][ivx]);
						
// 						// initialize current node
// 						node_t here = { 0 };
// 						here.k = k;
// 						here.x = x; here.y = y;
// 						here.vx = vx;
// 						here.ix = ((k == 0) ? 0 : ix); here.iy = ((k == 0) ? 0 : iy); 
// 						here.ivx = ((k == 0) ? 0 : ivx);
						
// 						if (here.ix > this->NX[k] || here.iy > this->NY[k] || here.ivx > this->NVX[k]) {
// 							fprintf(stderr, "invalid state on the grid\n"); 
// 							exit(1);
// 						}

// 						// final time step costs are zero
// 						if (k == this->numsteps) {
// 							if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
// 								J[here.k][here.ix][here.iy][here.ivx] = DBL_MAX;
// 							else
// 								J[here.k][here.ix][here.iy][here.ivx] = 0.0;			

// 							Ux[here.k][here.ix][here.iy][here.ivx] = 0.0;
// 							Uy[here.k][here.ix][here.iy][here.ivx] = 0.0;
							
// 							continue;
// 						}

// 						node_t next = { 0 };
// 						double cost, phi;
// 						double *Jhere, *Jnext, *Uxhere, *Uyhere;
						
// 						Jhere = &J[here.k][here.ix][here.iy][here.ivx];
// 						Uxhere = &Ux[here.k][here.ix][here.iy][here.ivx];
// 						Uyhere = &Uy[here.k][here.ix][here.iy][here.ivx];
						
// 						*Jhere = DBL_MAX;
// 						next.k = here.k + 1;

// 						// if (crash)
// 						if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
// 							continue;
						
// 						// for each control variable
// 						for (auto ux : UX) {
// 							for (auto uy : UY) {
// 								next.x = here.x + here.vx*this->step + 0.5*ux*pow(this->step, 2);
// 								next.y = here.y + uy*this->step;
// 								next.vx = here.vx + ux*this->step;

// 								// if ( feasible discrete point )
// 								if (discretizexv(&next, next.k) != 0)
// 									continue;

// 								// TODO: Interpolation
								
// 								// compute cost
// 								Jnext = &J[next.k][next.ix][next.iy][next.ivx];

// 								phi = 10.0 * 0.5*pow(ux, 2) + pow(here.vx - v.getDersiredSpeedX(), 2) + pow(uy - v.getDersiredSpeedY(), 2);
// 								cost =  phi + *Jnext;

// 								if (cost < *Jhere) {
// 									*Jhere = cost;
// 									*Uxhere = ux;
// 									*Uyhere = uy;
// 								}
// 							}					
// 						}
// 					}
// 				}
// 			}
// 		}
// 		this->printSolution(v, data, r, it);
// 		this->plotSolution(v, it);
// 		it_toc = clock();
// 		it_cpu_time = (double)(it_toc - it_tic) / CLOCKS_PER_SEC;
// 		fprintf(stderr, "Iteration %02d | Δux=Δvx: %.3f | Δx: %.4f | NX: %d | NY: %02d | NVX: %02d | NUX: %02d | NUY: %02d | J_cur: %04.4f | J_prev: %.4f | CPU time: %.3f \n", 
// 			it + 1, dux, dx, maxNX, maxNY, maxNVX, NUX, NUY, opt_J_it, opt_J_prev, it_cpu_time);

// 		// if (it >= 1 || this->dux0 < 0.125) { 
// 		if (this->dux0 <= 0.125 && this->duy0 <= 0.125) {
// 		// if (fabs(opt_J_it - opt_J_prev) < 0.1) {
// 			break;
// 		}
// 		if (opt_J_it == opt_J_prev) {
// 			this->dux0 /= 2.0; 
// 			this->duy0 /= 2.0;

// 			this->cx0 /= 2.0; this->cvx0 /= 2.0; this->cy0 /= 2.0;
// 		}
// 		opt_J_prev = opt_J_it;
		
// 		// re-adjust the initial trajectories
// 		for (int k = 0; k < numsteps + 1; k++) {
// 			v.init_x[k] = opt_x[k];
// 			v.init_y[k] = opt_y[k];
// 			v.init_vx[k] = opt_vx[k];
// 		}
// 		// clear memory and whatever else
// 		myFree();

// 		it++;
// 	}
// }

// void Controller::ddp(Vehicle v, Infrastructure r, Data data) {
// 	int it = 0;
// 	double x, y, vx;
// 	double x0 = v.getX0(), y0 = v.getY0(), vx0 = v.getVX0();
// 	v.setDersiredSpeedX(data.vdx); v.setDersiredSpeedY(data.vdy);

// 	clock_t it_tic, it_toc;
// 	double it_cpu_time;

// 	while(true) {
// 		it_tic = clock();
// 		// TODO: check this
// 		this->feasible_domain(v, r, it);
// 		this->allocations();
// 		// fprintf(stderr, "(%d) dux: %.2f -- duy: %.2f\n", it+1, this->dux, this->duy);

// 		for (int k = this->numsteps; k >= 0; k--) {
// 			// for each state variable
// 			for (int ix = 0; ix < ((k == 0) ? 1 : this->NX[k]); ix++) {
// 				x = ((k==0) ? x0 : this->X[k][ix]);
// 				for (int iy = 0; iy < ((k == 0) ? 1 : this->NY[k]); iy++) {
// 					y = ((k==0) ? y0 : this->Y[k][iy]);
// 					for (int ivx = 0; ivx < ((k == 0) ? 1 : this->NVX[k]); ivx++) {
// 						vx = ((k==0) ? vx0 : this->VX[k][ivx]);
						
// 						// initialize current node
// 						node_t here = { 0 };
// 						here.k = k;
// 						here.x = x; here.y = y;
// 						here.vx = vx;
// 						here.ix = ((k == 0) ? 0 : ix); here.iy = ((k == 0) ? 0 : iy); 
// 						here.ivx = ((k == 0) ? 0 : ivx);
						
// 						if (here.ix > this->NX[k] || here.iy > this->NY[k] || here.ivx > this->NVX[k]) {
// 							fprintf(stderr, "invalid state on the grid\n"); 
// 							exit(1);
// 						}

// 						// final time step costs are zero
// 						if (k == this->numsteps) {
// 							if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
// 								J[here.k][here.ix][here.iy][here.ivx] = DBL_MAX;
// 							else
// 								J[here.k][here.ix][here.iy][here.ivx] = 0.0;			

// 							Ux[here.k][here.ix][here.iy][here.ivx] = 0.0;
// 							Uy[here.k][here.ix][here.iy][here.ivx] = 0.0;
							
// 							continue;
// 						}

// 						node_t next = { 0 };
// 						double cost, phi;
// 						double *Jhere, *Jnext, *Uxhere, *Uyhere;
						
// 						Jhere = &J[here.k][here.ix][here.iy][here.ivx];
// 						Uxhere = &Ux[here.k][here.ix][here.iy][here.ivx];
// 						Uyhere = &Uy[here.k][here.ix][here.iy][here.ivx];
						
// 						*Jhere = DBL_MAX;
// 						next.k = here.k + 1;

// 						// if (crash)
// 						if (v.collision(data, here.x, here.y, here.vx, here.k) == 1)
// 							continue;
						
// 						// for each control variable
// 						for (auto ux : UX) {
// 							for (auto uy : UY) {
// 								next.x = here.x + here.vx*this->step + 0.5*ux*pow(this->step, 2);
// 								next.y = here.y + uy*this->step;
// 								next.vx = here.vx + ux*this->step;

// 								// if ( feasible discrete point )
// 								if (discretizexv(&next, next.k) != 0)
// 									continue;

// 								// TODO: Interpolation
								
// 								// compute cost
// 								Jnext = &J[next.k][next.ix][next.iy][next.ivx];

// 								phi = 10.0 * 0.5*pow(ux, 2) + pow(here.vx - v.getDersiredSpeedX(), 2) + pow(uy - v.getDersiredSpeedY(), 2);
// 								cost =  phi + *Jnext;

// 								if (cost < *Jhere) {
// 									*Jhere = cost;
// 									*Uxhere = ux;
// 									*Uyhere = uy;
// 								}
// 							}					
// 						}
// 					}
// 				}
// 			}
// 		}
// 		this->printSolution(v, data, r, it);
// 		this->plotSolution(v, it);
// 		it_toc = clock();
// 		it_cpu_time = (double)(it_toc - it_tic) / CLOCKS_PER_SEC;
// 		fprintf(stderr, "Iteration %02d | Δux=Δvx: %.3f | Δx: %.4f | NX: %d | NY: %02d | NVX: %02d | NUX: %02d | NUY: %02d | J_cur: %.4f | J_prev: %.4f | CPU time: %.3f \n", 
// 			it + 1, dux, dx, maxNX, maxNY, maxNVX, NUX, NUY, opt_J_it, opt_J_prev, it_cpu_time);

// 		// if (it >= 1 || this->dux0 < 0.125) { 
// 		if (this->dux0 <= 0.125 && this->duy0 <= 0.125) {
// 		// if (fabs(opt_J_it - opt_J_prev) < 0.1) {
// 			break;
// 		}
// 		if (opt_J_it == opt_J_prev) {
// 			this->dux0 /= 2.0; 
// 			this->duy0 /= 2.0;

// 			this->cx0 /= 2.0; this->cvx0 /= 2.0; this->cy0 /= 2.0;
// 		}
// 		opt_J_prev = opt_J_it;
		
// 		// re-adjust the initial trajectories
// 		for (int k = 0; k < numsteps + 1; k++) {
// 			v.init_x[k] = opt_x[k];
// 			v.init_y[k] = opt_y[k];
// 			v.init_vx[k] = opt_vx[k];
// 		}
// 		// clear memory and whatever else
// 		myFree();

// 		it++;
// 	}
// }

// void Controller::plotDomain(Vehicle v, int it) {

// 	if (it == 0) {
// 		string outPath = "outputs/domain";
// 		deleteDirectoryContents(outPath);
// 	}

// 	char fileBuffer[32];
//     snprintf(fileBuffer, sizeof(char) * 32, "outputs/domain/domain_%02i.py", it);
	
// 	FILE *fout;
// 	fout = fopen(fileBuffer, "w");

// 	/* python file */
// 	fputs("domain = {\n", fout);
	
// 	// initial x
// 	fputs("'x0': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", v.init_x[k]);
// 	fputs("],\n", fout);

// 	// initial y
// 	fputs("'y0': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", v.init_y[k]);
// 	fputs("],\n", fout);

// 	// initial vx
// 	fputs("'vx0': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", v.init_vx[k]);
// 	fputs("],\n", fout);
	
// 	// vx
// 	fputs("'vxmin': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", vxmin[k]);
// 	fputs("],\n", fout);

// 	fputs("'vxmax': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", vxmax[k]);
// 	fputs("],\n", fout);

// 	fputs("'vx': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) {
// 		fputs("\t[", fout);
// 		for (int ivx = 0; ivx < NVX[k]; ivx++)
// 			fprintf(fout, "%.4f,", VX[k][ivx]);
// 		fputs("],\n", fout);
// 	}
// 	fputs("],\n", fout);

// 	// x
// 	fputs("'xmin': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", xmin[k]);
// 	fputs("],\n", fout);

// 	fputs("'xmax': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", xmax[k]);
// 	fputs("],\n", fout);

// 	fputs("'x': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) {
// 		fputs("\t[", fout);
// 		for (int ivx = 0; ivx < NX[k]; ivx++)
// 			fprintf(fout, "%.4f,", X[k][ivx]);
// 		fputs("],\n", fout);
// 	}
// 	fputs("],\n", fout);

// 	// y
// 	fputs("'ymin': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", ymin[k]);
// 	fputs("],\n", fout);

// 	fputs("'ymax': [", fout);
// 	for (int k = 0; k < numsteps; k++)
// 		fprintf(fout, "%.4f,", ymax[k]);
// 	fputs("],\n", fout);

// 	fputs("'y': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) {
// 		fputs("\t[", fout);
// 		for (int ivx = 0; ivx < NY[k]; ivx++)
// 			fprintf(fout, "%.4f,", Y[k][ivx]);
// 		fputs("],\n", fout);
// 	}
// 	fputs("],\n", fout);

// 	fprintf(fout, "'k':%d,\n", numsteps);
// 	fprintf(fout, "'Step':%f,\n", step);
// 	fprintf(fout, "'iter':%d,\n", it + 1);
// 	fprintf(fout, "'dux':%f,\n", dux);
// 	fprintf(fout, "'duy':%f,\n", duy);
// 	fprintf(fout, "'dvx':%f,\n", dvx);
// 	fprintf(fout, "'dx':%f,\n", dx);
// 	fprintf(fout, "'dy':%f,\n", dy);
// 	fprintf(fout, "}\n\n");


// 	fclose(fout);
// }

// void Controller::plotSolution(Vehicle v, int it) {

// 	if (it == 0) {
// 		string outPath = "outputs/solution";
// 		deleteDirectoryContents(outPath);
// 	}

// 	char fileBuffer[32];
//     snprintf(fileBuffer, sizeof(char) * 32, "outputs/solution/solution_%02i.py", it);
	
// 	FILE *fout;
// 	fout = fopen(fileBuffer, "w");

// 	/* python file */
// 	fputs("solution = {\n", fout);
	
// 	// vx
// 	fputs("'vx': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) 
// 		fprintf(fout, "%.4f,", opt_vx[k]);
// 	fputs("],\n", fout);

// 	// x
// 	fputs("'x': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) 
// 		fprintf(fout, "%.4f,", opt_x[k]);
// 	fputs("],\n", fout);

// 	// y
// 	fputs("'y': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) 
// 		fprintf(fout, "%.4f,", opt_y[k]);
// 	fputs("],\n", fout);

// 	// ux
// 	fputs("'ux': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) 
// 		fprintf(fout, "%.4f,", opt_ux[k]);
// 	fputs("],\n", fout);

// 	// uy
// 	fputs("'uy': [\n", fout);
// 	for (int k = 0; k < numsteps; k++) 
// 		fprintf(fout, "%.4f,", opt_uy[k]);
// 	fputs("],\n", fout);

// 	fprintf(fout, "'k':%d,\n", numsteps);
// 	fprintf(fout, "'Step':%f,\n", step);
// 	fprintf(fout, "'iter':%d,\n", it + 1);
// 	fprintf(fout, "}\n\n");


// 	fclose(fout);
// }
