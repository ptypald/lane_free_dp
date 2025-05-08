#include "dddp.h"
int UX_CORRIDOR = 1;

void DDDP::initial_trajectory() {
    if (init_x.size() != numsteps)
        fprintf(stderr, "re-adjust initial trajectories: %ld --> %d\n", init_x.size(), numsteps + 1);
    
    int adjust = (int)(init_x.size() / (numsteps + 1)) + 1;

    for (int i = 0; i < init_x.size(); i++) {
        if (i % adjust == 0)
            continue;

        init_x[i] = DBL_MAX;
        init_vx[i] = DBL_MAX;
        init_y[i] = DBL_MAX;
    }

    init_x.erase(remove(init_x.begin(), init_x.end(), DBL_MAX), init_x.end());
    init_vx.erase(remove(init_vx.begin(), init_vx.end(), DBL_MAX), init_vx.end());
    init_y.erase(remove(init_y.begin(), init_y.end(), DBL_MAX), init_y.end());
    
	// create intitial control trajectories
	for (int i = 1; i < numsteps + 1; i++) {
		init_ux.push_back((init_vx[i] - init_vx[i - 1]) / step);
		init_uy.push_back((init_y[i] - init_y[i - 1]) / step);
	}



	// print trajectories for mathematica
	// for (int i = 1; i < numsteps; i++)
	// 	fprintf(stderr, "%.7f, ", (init_uy[i] - init_uy[i - 1]) / step);
	// fprintf(stderr, "\n");

}

void DDDP::initialize_parameters(Vehicle veh) {
    dux0 = 1.0;
    uxmin = veh.UXMIN; uxmax = veh.UXMAX;
    duy0 = 1.0;
    uymin = veh.UYMIN; uymax = veh.UYMAX;

    vxMin = veh.VXMIN; vxMax = veh.VXMAX;
    yMin = 0.0; yMax = 10.0;
    xMin = x0; xMax = DBL_MAX;

    // cvx0 = (2.0)/step; cx0 = (2.0*cvx0)/step; cy0 = (2.0*duy0)/step;	// how many point around the initial trajectories (number indicates point only on one side, so totall points are x2 + 1)
	cvx0 = (2.0)/step; cx0 = (2.0*cvx0)/step; cy0 = (2.0*duy0)/step;
}

void DDDP::feasible_domain(int it) {
	vector<double> temp, temp2;
	double mod_x, mod_y;

	// clear vectors
	this->UX.clear(); this->UY.clear();
	this->vxmin.clear(); this->vxmax.clear(); this->NVX.clear(); this->VX.clear();
	this->xmin.clear(); this->xmax.clear(); this->NX.clear(); this->X.clear();
	this->ymin.clear(); this->ymax.clear(); this->NY.clear(); this->Y.clear();
	UXK.clear(); UYK.clear();

	/* ux and uy domain */
	this->dux = this->dux0; this->duy = this->duy0;

	// TODO: add corridor for the acceleration trajectories
	if (UX_CORRIDOR == 0) {
		this->NUX = (int)ceil((this->uxmax - this->uxmin) / this->dux) + 1;
		this->NUY = (int)ceil((this->uymax - this->uymin) / this->duy) + 1;
		for (int k = 0; k < this->numsteps + 1; k++) {
			for (int iux = 0; iux < this->NUX; iux++) {
				temp.push_back(this->uxmin + this->dux * iux);
			}
			for (int iuy = 0; iuy < this->NUY; iuy++) {
				temp2.push_back(this->uymin + this->duy * iuy);
			}
			UXK.push_back(temp); UYK.push_back(temp2);
		}
		temp.clear(); temp2.clear();
	}
	else {
		this->NUX = (int)((2.0*dux) / this->dux) + 1;
		this->NUY = (int)((2.0*duy) / this->duy) + 1;
		for (int k = 0; k < this->numsteps; k++) {
			mod_x = fmod(init_ux[k], dux);
			mod_y = fmod(init_uy[k], dux);
			
			if (mod_x != 0.0 && it == 0) {
				if (init_ux[k] < 0.0)
					init_ux[k] = MAX(init_ux[k] - (1.0 + mod_x), uxmin);
				else
					init_ux[k] = MIN(init_ux[k] + (1.0 - mod_x), uxmax);
			}
			if (mod_y != 0.0 && it == 0) {
				if (init_uy[k] < 0.0)
					init_uy[k] = MAX(init_uy[k] - (1.0 + mod_y), uymin);
				else
					init_uy[k] = MIN(init_uy[k] + (1.0 - mod_y), uymax);
			}

			double test, test2;
			for (int iux = 0; iux < this->NUX; iux++) {
				test = MAX(((init_ux[k]) - dux) + (dux * iux), uxmin);
				test2 = MIN(test, uxmax);
				temp.push_back(test2);

				test = MAX(((init_uy[k]) - duy) + (duy * iux), uymin);
				test2 = MIN(test, uymax);
				temp2.push_back( test2 );
			}
			
			UXK.push_back(temp); UYK.push_back(temp2);
			temp.clear(); temp2.clear();

			sort( UXK[k].begin(), UXK[k].end() );
			UXK[k].erase( unique( UXK[k].begin(), UXK[k].end() ), UXK[k].end() );

			sort( UYK[k].begin(), UYK[k].end() );
			UYK[k].erase( unique( UYK[k].begin(), UYK[k].end() ), UYK[k].end() );
		}
		
	}

	/* vx domain */
	this->dvx = this->dux * this->step; 
	this->cvx = this->cvx0 * this->dvx;
	for (int k = 0; k < this->numsteps + 1; k++) {
		mod_x = fmod(init_vx[k], dvx);
		if (it == 0 && mod_x != 0.0) {
			init_vx[k] = init_vx[k] + (1.0 - mod_x);
			this->vxmin.push_back(MAX(init_vx[k] - this->cvx, vxMin));
			this->vxmax.push_back(MIN(init_vx[k] + this->cvx, vxMax));
		}
		else {
			this->vxmin.push_back(MAX(init_vx[k] - this->cvx, vxMin));
			this->vxmax.push_back(MIN(init_vx[k] + this->cvx, vxMax));
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
	cx = cx0 * dx;
	for (int k = 0; k < this->numsteps + 1; k++) {
		mod_x = fmod(init_x[k], dx);
		if (it == 0 && mod_x != 0.0) {
			init_x[k] = init_x[k] + (1.0 - mod_x);
			this->xmin.push_back(MAX(init_x[k] - this->cx, xMin));
			this->xmax.push_back(MIN(init_x[k] + this->cx, xMax));
		}
		else {
			this->xmin.push_back(MAX(init_x[k] - this->cx, xMin));
			this->xmax.push_back(MIN(init_x[k] + this->cx, xMax));
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
	this->dy = this->duy * this->step;
	this->cy = this->cy0 * this->dy;
	for (int k = 0; k < this->numsteps + 1; k++) {
		mod_x = fmod(init_y[k], dy);
		if (it == 0 && mod_x != 0.0) {
			init_y[k] = init_y[k] + (1.0 - mod_x);
			this->ymin.push_back(MAX(init_y[k] - this->cy, yMin + 1.0)); // +1 and -1 for the dimensions of the vehicle (TODO: add the real dimensions)
			this->ymax.push_back(MIN(init_y[k] + this->cy, yMax - 1.0));
		}
		else {
			this->ymin.push_back(MAX(init_y[k] - this->cy, yMin + 1.0));
			this->ymax.push_back(MIN(init_y[k] + this->cy, yMax - 1.0));
		}
		this->NY.push_back((int)ceil((this->ymax[k] - this->ymin[k]) / this->dy) + 1);

		for (int iy = 0; iy < this->NY[k]; iy++)
			temp.push_back(this->ymin[k] + this->dy * iy);
		this->Y.push_back(temp);
		temp.clear();
	}
	maxNY = *max_element(NY.begin(), NY.end());
	
	plotDomain(it);
}

int DDDP::discretizexv(node_t *state, int k) {
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

int DDDP::collision(Vehicle veh, double x, double y, double vx, int k){
    double obsx, obsy, obsvx;
    for (int i = 0; i < numObs; i++) {
        obsx = obs[i].x[k], obsy = obs[i].y[k], obsvx = obs[i].vx[k];
        if (veh.collision_per_obs(x, obsx, y, obsy, vx, obsvx) == 1)
             return 1;                  
    }
    return 0;
}

void DDDP::run(Vehicle v) {
	int it = 0;
	double x, y, vx;
    int collided;

	clock_t it_tic, it_toc, all_tic, all_toc;
	double it_cpu_time, total_cpu_time;

	all_tic = clock();
	while(true) {
		it_tic = clock();
		
		feasible_domain(it); allocations();
		for (int k = numsteps; k >= 0; k--) {
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
						if (here.k == numsteps) {                           
                            collided = 0;
                            
                            // for (int o = 0; o < numObs; o++){
							if (collision(v, here.x, here.y, here.vx, here.k) == 1) {
								collided = 1;
								// break;
							}
                                
                            // }

							if (collided)
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

						if (collision(v, here.x, here.y, here.vx, here.k) == 1)
							continue;
					 
						// for each control variable
						for (auto ux : UXK[here.k]) {
							for (auto uy : UYK[here.k]) {
								next.x = here.x + here.vx*this->step + 0.5*ux*pow(this->step, 2);
								next.y = here.y + uy*this->step;
								next.vx = here.vx + ux*this->step;

								// if ( feasible discrete point )
								if (discretizexv(&next, next.k) != 0)
									continue;

								// TODO: Interpolation
								
								// compute cost
								Jnext = &J[next.k][next.ix][next.iy][next.ivx];

								phi = 10.0 * 0.5*pow(ux, 2) + pow(here.vx - vdx, 2) + pow(uy - 0.0, 2);
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
		this->forwardRun(it);
		this->plotSolution(it);

		it_toc = clock();
		it_cpu_time = (double)(it_toc - it_tic) / CLOCKS_PER_SEC;
		if (DISPLAY)
			fprintf(stderr, "Iteration %02d | Δux=Δvx: %.3f | Δx: %.4f | NX: %d | NY: %02d | NVX: %02d | NUX: %02d | NUY: %02d | J_cur: %04.4f | J_prev: %.4f | CPU time: %.3f \n", 
				it + 1, dux, dx, maxNX, maxNY, maxNVX, NUX, NUY, opt_J_it, opt_J_prev, it_cpu_time);
 
		if (opt_J_it == opt_J_prev) {
			this->dux0 /= 2.0; 
			this->duy0 /= 2.0;
		}
		opt_J_prev = opt_J_it;

		if ( (this->dux0 < 0.125 && this->duy0 < 0.125) ) {
		// if (fabs(opt_J_it - opt_J_prev) < 0.1) {
			break;
		}
		
		// re-adjust the initial trajectories
		for (int k = 0; k < numsteps + 1; k++) {
			init_x[k] = opt_x[k];
			init_y[k] = opt_y[k];
			init_vx[k] = opt_vx[k];
	
			init_ux[k] = opt_ux[k];
			init_uy[k] = opt_uy[k];
		}
   
		// clear memory and whatever else
		myFree();

		it++;
	}
	all_toc = clock();
	total_cpu_time = (double)(all_toc - all_tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- (DDDP) -- Total Iterations %02d | J_opt: %04.4f | Total CPU time: %.3f seconds \n", 
		it + 1, opt_J_it, total_cpu_time);
}

void DDDP::forwardRun(int it){

	node_t *nextptr;
	node_t *state;
	double ux, uy, real_cost = 0.0;

	opt_ux.clear(); opt_uy.clear();
	opt_x.clear(); opt_y.clear(); opt_vx.clear();

	node_t initial = { 0 };
	initial.k = 0;
	initial.x = x0;
	initial.y = y0;
	initial.vx = vx0; //getVX0();
	initial.ix = initial.iy = initial.ivx = 0;

	state = (node_t*)calloc(sizeof(node_t), this->numsteps + 1);
	state[0] = initial;

	if (DISPLAY) {
		fprintf(stderr, "k \t x \t\t y \t\t vx \t\t ux \t\t uy \t\t J \t\t real cost \n");
		fprintf(stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	}
	for (int k = 0; k <= this->numsteps; k++){
		ux = ACCESS(Ux, state[k]); opt_ux.push_back(ux);
		uy = ACCESS(Uy, state[k]); opt_uy.push_back(uy);
		opt_x.push_back(state[k].x); opt_y.push_back(state[k].y); opt_vx.push_back(state[k].vx);

		// real cost
		real_cost += 0.5*pow(ux, 2) + 0.5*pow(state[k].vx - vdx, 2) + 0.5*pow(uy - 0.0, 2);
		
		if (DISPLAY) {
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
	exportSimData(it);

}

void DDDP::plotDomain(int it) {
	char fileBuffer[64];
    snprintf(fileBuffer, sizeof(char) * 64, "../outputs/dddp/domain/domain_%02i.py", it);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");
	/* python file */
	fputs("domain = {\n", fout);
	
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

	// ux
	fputs("'ux': [\n", fout);
	for (int k = 0; k < numsteps; k++) {
		fputs("\t[", fout);
		for (int iux = 0; iux < NUX; iux++)
			fprintf(fout, "%.4f,", UXK[k][iux]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);

	fputs("'uy': [\n", fout);
	for (int k = 0; k < numsteps; k++) {
		fputs("\t[", fout);
		for (int iuy = 0; iuy < NUY; iuy++)
			fprintf(fout, "%.4f,", UYK[k][iuy]);
		fputs("],\n", fout);
	}
	fputs("],\n", fout);
	// }

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

void DDDP::plotSolution(int it) {

	char fileBuffer[64];
    snprintf(fileBuffer, sizeof(char) * 64, "../outputs/dddp/solution/solution_%02i.py", it);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");

	/* python file */
	fputs("solution = {\n", fout);
	
	// initial x
	fputs("'x0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_x[k]);
	fputs("],\n", fout);

	// initial y
	fputs("'y0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_y[k]);
	fputs("],\n", fout);

	// initial vx
	fputs("'vx0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_vx[k]);
	fputs("],\n", fout);

	// initial ux
	fputs("'ux0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_ux[k]);
	fputs("],\n", fout);

	// initial ux
	fputs("'uy0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_uy[k]);
	fputs("],\n", fout);

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

void DDDP::allocations() {

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

void DDDP::myFree() {
	
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