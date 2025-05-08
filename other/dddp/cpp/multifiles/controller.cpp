#include "controller.h"
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.ivx]

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

void Controller::myAlloc() {

	J = (double ****)calloc(sizeof(double***), numsteps + 1);
	Ux = (double ****)calloc(sizeof(double***), numsteps + 1);
	Uy = (double ****)calloc(sizeof(double***), numsteps + 1);
	for (int k = 0; k < numsteps + 1; k++) {
		J[k] = (double ***)calloc(sizeof(double**), maxNX);
		Ux[k] = (double ***)calloc(sizeof(double**), maxNX);
		Uy[k] = (double ***)calloc(sizeof(double**), maxNX);
		for (int ix = 0; ix < maxNX; ix++) {
			J[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			Ux[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			Uy[k][ix] = (double **)calloc(sizeof(double*), maxNY);
			for (int iy = 0; iy < maxNY; iy++) {
				J[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
				Ux[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
				Uy[k][ix][iy] = (double *)calloc(sizeof(double), maxNVX);
			}
		}
	}
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

int Controller::collision(Vehicle v, double here_x, double here_y, double here_vx, int here_k) {

	for (int i = 0; i < v.obs_n; i++) {
		if ( (here_x < v.obs[i].x[here_k] + here_vx*v.safety_x && here_x > v.obs[i].x[here_k] - v.obs[i].vx[here_k]*v.safety_x) && (here_y < v.obs[i].y[here_k] + 2.0 && here_y > v.obs[i].y[here_k] - 2.0) ) {
			return 1;
		}
	}
	
	return 0;
} 

void Controller::printSolution(Vehicle v, int it){

	node_t *nextptr;
	node_t *state;
	double ux, uy, real_cost = 0.0;

	opt_x.clear(); opt_y.clear(); opt_vx.clear();

	node_t initial = { 0 };
	initial.k = 0;
	initial.x = v.getX0();
	initial.y = v.getY0();
	initial.vx = v.getVX0();
	initial.ix = initial.iy = initial.ivx = 0;

	state = (node_t*)calloc(sizeof(node_t), this->numsteps + 1);
	state[0] = initial;

	// if (data.DISPLAY) {
		fprintf(stderr, "k \t x \t\t y \t\t vx \t\t ux \t\t uy \t\t J \t\t real cost \n");
		fprintf(stderr, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	// }
	for (int k = 0; k <= this->numsteps; k++){
		ux = ACCESS(Ux, state[k]); uy = ACCESS(Uy, state[k]);
		opt_x.push_back(state[k].x); opt_y.push_back(state[k].y); opt_vx.push_back(state[k].vx);

		// real cost
		real_cost += 0.5*pow(ux, 2) + 0.5*pow(state[k].vx - v.getDersiredSpeedX(), 2) + 0.5*pow(uy - v.getDersiredSpeedY(), 2);
		
		// if (data.DISPLAY) {
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
		// }

		if (k < numsteps) {
			if ((nextptr = getnext(state[k], ux, uy)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}
	}

	opt_J_it = ACCESS(J, state[0]);
	free(state);
	// data.exportSimData(v, inf, it);

}



