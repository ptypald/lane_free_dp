#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>

// for file creation
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#define Power(a, b) pow(a, b)
// #define fopen fopen_s
#define sscanf_s sscanf
//#define fopen_s fopen
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iv]

#define STOCHASTIC 1
#define TRIANGULAR 0
#define FIXED_TIME 0

#define TIME_STEP 1.0

#define DELTA_V 4.0
#define DELTA_X_mult 5.0

#define TERM_CRITERION 0.1
#define TERM_CRITERION_2 0.00001
#define TERM_CRITERION_DU 1

#define ALLOC_STATES 124

/* modes */
#define DISPLAY 0
#define OUTPUT 1
#define OBSTACLES 0

double du, init_du = 0.5;

double XMAX = 150.0, XMIN = 0.0;
double VMAX = 16.0,  VMIN = 0.0;
double T_min = 10.0, T_max = 30.0, T_tl = 0.0;  	// possible change times from T_min to T_max
double *x_min, *x_max;								// x domain (difDP)
double u_min = -3.0, u_max = 3.0;					// u domain

typedef struct node {
	int k;
	int ix, iv;
	double x, v;
	double _x, _v;
} node_t;

typedef struct {
	double u;
} ctrl_t;

double ***J, ***J_opt, ***escape, ***te_opt;
ctrl_t ***C;
double *p;

/* new */
double *opt_x, *opt_v, *opt_u;
double opt_J_it;
double *init_x, *init_v, *init_u;
double *initial_x, *initial_v, *initial_u;
// double *initial_x, *initial_v, *initial_a; //fix all the vectors, check for un-needed
/* UP problem */
double te_up;
double Dv, Dx;

/* problem */
struct {
	double T;
	double x0, y0, v0;
	int numsteps;
	double ***opt_te;

	/* obstacles */
	double safety;
	int n;
}P = { 0 };

/* state variable domains */
struct {
	int *NX, *NXlive_min, *NXlive_max;
	double dx, **X;
	double *xmin, *xmax;

	int *NV, *NVlive_min, *NVlive_max;
	double *vmin, *vmax, dv, **V;

	int NU;
	double umin, umax, du;
}D = { 0 };

/* obstacles */
struct {
	double **obst_x;
	double **obst_y;
	double **obst_v;
}O = { 0 };

static node_t * getnext(node_t here, ctrl_t c);

static int discretizexv(node_t *state, int k);

static double crash(node_t state) {
	int i;
	double dx;
	double sgap;
	double cost = 0;
	double L = 5.0;

	for (i = 0; i < P.n; i++) {
		dx = fabs(O.obst_x[state.k][i] - state.x);

		if (O.obst_x[state.k][i] - state.x >= 0) {
			sgap = (state.v * (P.safety) + L);
		}
		else {
			sgap = (O.obst_v[state.k][i] * (P.safety) + L);
		}
		if (dx <= sgap){
			cost += 1.0;
		}

	}
	return cost;
}

static void init_S(int numsteps, double T, int n, double obst_x, double obst_v)
{
	int i, k;

	/* allocate */
	assert((O.obst_x = (double**)calloc(sizeof(double*), numsteps + 1)));
	assert((O.obst_v = (double**)calloc(sizeof(double*), numsteps + 1)));
	for (k = 0; k < numsteps + 1; k++) {
		assert((O.obst_x[k] = (double*)calloc(sizeof(double), n)));
		assert((O.obst_v[k] = (double*)calloc(sizeof(double), n)));
	}

	for (i = 0; i < n; i++) {
		O.obst_x[0][i] = obst_x;
		O.obst_v[0][i] = obst_v;
	}

	/* compute */
	for (k = 0; k < P.numsteps; k++) {
		for (i = 0; i < n; i++) {
			O.obst_x[k + 1][i] = O.obst_x[k][i] + O.obst_v[k][i] * P.T;
			O.obst_v[k + 1][i] = O.obst_v[k][i];
		}
	}
}

static void init_D(double T, int iteration) {

	int ix, iv;
	int i, k;
  
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	assert((D.V = (double**)calloc(sizeof(double*), ALLOC_STATES)));
	for (i = 0; i < ALLOC_STATES; i++)
		assert((D.V[i] = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.X = (double**)calloc(sizeof(double*), ALLOC_STATES)));
	for (i = 0; i < ALLOC_STATES; i++)
		assert((D.X[i] = (double*)calloc(sizeof(double), P.numsteps)));
  
	assert((D.vmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.vmin = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.xmax = (double*)calloc(sizeof(double), P.numsteps)));
	assert((D.xmin = (double*)calloc(sizeof(double), P.numsteps)));

	assert((D.NV = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NX = (int*)calloc(sizeof(int), P.numsteps)));
  
	/* u domain: du = dv/T */
	D.umin = u_min;
	D.umax = u_max;

	D.du = du;
	D.NU = (int)ceil((D.umax - D.umin) / D.du) + 1;
  
	/* v domain */
	D.dv = D.du;
	Dv = DELTA_V;

	for (k = 0; k < P.numsteps; k++){		
        D.vmin[k] = MAX(init_v[k] - Dv * D.dv, VMIN);
        D.vmax[k] = MIN(init_v[k] + Dv * D.dv, VMAX);
	
		D.NV[k] = (int)ceil((D.vmax[k] - D.vmin[k]) / D.dv) + 1;
		
		for (iv = 0; iv < D.NV[k]; iv++) {
			D.V[iv][k] = D.vmin[k] + D.dv * iv;
		// fprintf(stderr, "D.V[%d][%d]: %.3f \n", iv, k, D.V[iv][k]);
		}
	}
	
	/* x domain */
	D.dx = 0.5 * D.du;
	Dx = DELTA_V * DELTA_X_mult;
	for (k = 0; k < P.numsteps; k++) {
		if (iteration == 0) {
			D.xmin[k] = floor(MAX(init_x[k] - Dx * D.dx, XMIN));
			D.xmax[k] = ceil(MIN(init_x[k] + Dx * D.dx, XMAX));
	
		}
		else {
			D.xmin[k] = MAX(init_x[k] - Dx * D.dx, XMIN);
			D.xmax[k] = MIN(init_x[k] + Dx * D.dx, XMAX);

		}
		D.NX[k] = (int)ceil((D.xmax[k] - D.xmin[k]) / D.dx) + 1;

		for (ix = 0; ix < D.NX[k]; ix++) {
			D.X[ix][k] = D.xmin[k] + D.dx * ix;
		// fprintf(stderr, "(%d) D.X[%d][%d]: %.3f \n", D.NX[k], ix, k, D.X[ix][k]);
		}
	}
  
	assert((D.NXlive_min = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NXlive_max = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NVlive_min = (int*)calloc(sizeof(int), P.numsteps)));
	assert((D.NVlive_max = (int*)calloc(sizeof(int), P.numsteps)));
	
	for (k = 0; k < P.numsteps; k++) {
		D.NXlive_min[k] = 0;
		D.NXlive_max[k] = D.NX[k];

		D.NVlive_min[k] = 0;
		D.NVlive_max[k] = D.NV[k];
	}
}

static void init_JC(void) {

	assert((J = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((C = (ctrl_t ***)calloc(sizeof(ctrl_t*), P.numsteps)));
	assert((J_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((te_opt = (double ***)calloc(sizeof(double*), P.numsteps)));
	assert((escape = (double ***)calloc(sizeof(double*), P.numsteps)));
	int k;
	for (k = 0; k < P.numsteps; k++) {
		assert((J[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((C[k] = (ctrl_t **)calloc(sizeof(ctrl_t*), ALLOC_STATES)));
		assert((J_opt[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((escape[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		assert((te_opt[k] = (double **)calloc(sizeof(double*), ALLOC_STATES)));
		int ix;
		for (ix = 0; ix < ALLOC_STATES; ix++) {
			assert((J[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((C[k][ix] = (ctrl_t *)calloc(sizeof(ctrl_t), ALLOC_STATES)));
			assert((J_opt[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((escape[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
			assert((te_opt[k][ix] = (double *)calloc(sizeof(double), ALLOC_STATES)));
		}
	}
	assert((p = (double *)calloc(sizeof(double), P.numsteps)));
}

static void free_JC(void) {
	int k, ix;
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < ALLOC_STATES; ix++) {
			free(J[k][ix]);
			free(C[k][ix]);
			free(J_opt[k][ix]);
			free(escape[k][ix]);
			free(te_opt[k][ix]);
		}
	}

	for (k = 0; k < P.numsteps; k++) {
		free(J[k]);
		free(C[k]);
		free(J_opt[k]);
		free(escape[k]);
		free(te_opt[k]);
	}
	free(J);
	free(C);
	free(J_opt);
	free(escape);
	free(te_opt);
}

static int discretizexv(node_t *state, int k) {
	int ix, iv;

	ix = (int)round((state->x - D.X[0][k]) / D.dx);
	iv = (int)round((state->v - D.V[0][k]) / D.dv);
	//fprintf(stderr, "ix: %d -- x: %.3f -- D.X[0][%d]: %.3f \n", ix, state->x, k, D.X[0][k]);
	if (!(ix >= 0 && ix < D.NX[k])) 
		return 1;
	
	if (!(iv >= 0 && iv < D.NV[k]))  
		return 2;

	state->x = D.X[(state->ix = ix)][k];
	state->v = D.V[(state->iv = iv)][k];
	//fprintf(stderr, "next.k: %d -- next.x: %.3f -- next.v: %.3f -- D.X[%d][%d]: %.3f \n", k, state->x, state->v, ix, k, D.X[ix][k]);
	return 0;
}

static node_t * getnext(node_t here, ctrl_t c) {
	static node_t next;
	int err;

	next.k = here.k + 1;
	next.x = next._x = here.x + here.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = next._v = here.v + c.u*P.T;
	if ((err = discretizexv(&next, next.k))) 
		return NULL;

	return &next;
}

static void printsol(node_t initial, int it) {
	
	node_t *nextptr;
	int k = 0;
	node_t *state;

	assert((state = (node_t*)calloc(sizeof(node_t), P.numsteps)));
	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c;
	
	FILE *f;
	if (DISPLAY) {
		fprintf(stderr, "\nk \t Position \t Speed    \t Control  \t Cost \t \t Hybrid  \t J_opt  \t te \t Probability \t Real Cost (1/2u^2) \t Real+escape \t Escape \n");
		fprintf(stderr, "-----------------------------------------------------------------------------------------------------------------\n");
	}
	if (OUTPUT) {
		double deltav = DELTA_V;
		char buf_dir[0x200];
		snprintf(buf_dir, sizeof(buf_dir), "x0=%.1f-v0=%.1f/", P.x0, P.v0);
		if (mkdir(buf_dir, 0777) == -1) {
			// fprintf(stderr, "Error - directory not created \n");
			// exit(1);
		}
		else
        	fprintf(stderr, "directory created \n");

		char buf[0x200];
		snprintf(buf, sizeof(buf), "(%02d) x0=%.1f - v0=%.1f - delta=%.1f - du=%.3f.txt", it + 1, P.x0, P.v0, DELTA_V, du);
		f = fopen(buf,"w");
		// fprintf(f, "k  Position  Speed  Control  Cost  Hybrid  J_opt  te  Probability  Real_Cost_(1/2u^2)  Real+escape   Escape \n");
	}
	
	double hybrid_cost = 0.0;
	double real_cost = 0.0;
	double real_escape_cost = 0.0; // 0.5u^2 + the escape cost in each stage
	for (k = 0; k < P.numsteps; k++) {
		c = ACCESS(C, state[k]);
		
		/* real cost*/
		real_cost += 0.5 * pow(c.u, 2);

		/* hybrid cost */
		if (k < T_min){
			hybrid_cost += 0.5 * pow(c.u, 2);
			real_escape_cost += 0.5 * pow(c.u, 2);
		}
		else{
			hybrid_cost += (1.0 - p[k]) * 0.5 * pow(c.u, 2) + p[k] * ACCESS(J_opt, state[k]);
			if (STOCHASTIC)
				real_escape_cost = real_cost + ACCESS(J_opt, state[k]);
			else {
				if (k == P.numsteps - 1)
					real_escape_cost = real_cost + ACCESS(J_opt, state[k]);	
				else
					real_escape_cost = real_cost + ACCESS(escape, state[k]);
			}
		}
		
		if (DISPLAY) {
			// fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n",
			fprintf(stderr, "%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n",
				k,
				state[k].x,
				state[k].v,
				c.u,
				ACCESS(J, state[k]),
				hybrid_cost,
				ACCESS(J_opt, state[k]),
				ACCESS(te_opt, state[k]),
				p[k],
				real_cost,
				real_escape_cost
				// ACCESS(escape, state[k]),
				// D.vmin[k],
				// D.vmax[k],
				// D.xmin[k],
				// D.xmax[k],
				// init_v[k],
				// init_x[k]
				);
		}
		
		if (OUTPUT) {
			fprintf(f, "%d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n",
				k,
				state[k].x,
				state[k].v,
				c.u,
				ACCESS(J, state[k]),
				hybrid_cost,
				ACCESS(J_opt, state[k]),
				ACCESS(te_opt, state[k]),
				p[k],
				real_cost,
				real_escape_cost,
				ACCESS(escape, state[k]),
				D.vmin[k],
				D.vmax[k],
				D.xmin[k],
				D.xmax[k],
				init_v[k],
				init_x[k],
				init_u[k]
				);		
		}
		
		if (k < P.numsteps - 1) {
			if ((nextptr = getnext(state[k], c)) == NULL) {
				fprintf(stderr, "invalid \n");
				exit(1);
			}
			state[k + 1] = *nextptr;
		}

		opt_x[k] = state[k].x;
		opt_v[k] = state[k].v;
		opt_u[k] = c.u;
	}
	
	opt_J_it = ACCESS(J, state[0]);
	free(state);
	if (OUTPUT)
		fclose(f);
}

static void anal_solution() {
  
  double t0 = T_max;
	double te;
	double x0 = opt_x[P.numsteps-1], v0 = opt_v[P.numsteps-1];
	
  fprintf(stderr, "x0: %.4f \t v0: %.4f \n", x0, v0);

  int ix = (int)round((x0 - D.X[0][P.numsteps-1]) / D.dx);
	int iv = (int)round((v0 - D.V[0][P.numsteps-1]) / D.dv);

  // te = ACCESS(te_opt, state[P.numsteps-1]);
  te = te_opt[P.numsteps-1][ix][iv];

  double xk = x0, vk = v0;
	int opNumsteps = (int)(te - t0) + 1;
	//double step = (te - t0) / (opNumsteps);
	int count = 0;
	double a_op[100];

  FILE *outOC;
  outOC = fopen("output_oc.txt", "w+");

  fprintf(stderr, "t0: %.4f \t te: %.4f \t opNum: %d \n", t0, te, opNumsteps);

	//for (double t = t0; t <= (te); t = t + step) {
	double step = (te - t0) / opNumsteps;
	for (int i = 0; i < opNumsteps; i++) {
		double t = (i*(te - t0) / (opNumsteps - 1) + t0);
		a_op[count] = UP_control(t, t0, x0, v0, te);
		fprintf(outOC, "%.4f \t", a_op[count]);
		count++;
	}
	
	for (int i = 0; i < count; i++) {
		xk += vk * step + 0.5*pow(step, 2)*a_op[i];
		vk += step * a_op[i];

		fprintf(outOC, "\n %.4f", xk);
		fprintf(outOC, "\t %.4f", vk);
	}

  fprintf(stderr, "hello!!\n");

}

static void dp(void) {
	
	int k, ix, iv;
	double x, v;
	double te;
	if (DISPLAY) {
		if (STOCHASTIC)
			fprintf(stderr, "Stochastic DP \n");
		else
			fprintf(stderr, "Deterministic DP \n");
		fprintf(stderr, "NU: %d -- NV: %d -- NX: %d \n", D.NU, D.NV[2], D.NX[2]);
	}
	//fprintf(stderr, "NU: %d -- NV: %d -- NX: %d \n", D.NU, D.NV[2], D.NX[2]);
	/* for each stage k */
	for (k = P.numsteps - 1; k >= 0; k--) {
		if (DISPLAY)
			fprintf(stderr, "%d ", k);
		/* for each possible x */
		for (ix = D.NXlive_min[k]; ix < ((k == 0) ? 1 : D.NXlive_max[k]); ix++) {
			x = ((k == 0) ? P.x0 : D.X[ix][k]);
			/* for each possible v */
			for (iv = D.NVlive_min[k]; iv < ((k == 0) ? 1 : D.NVlive_max[k]); iv++) {
				v = ((k == 0) ? P.v0 : D.V[iv][k]);

				node_t here = { 0 };
				here.k = k;
				here.x = x;
				here.v = v;
				here.ix = ((k == 0) ? 0 : ix);
				here.iv = ((k == 0) ? 0 : iv);
				
				if (here.ix > D.NX[k] || here.iv > D.NV[k]) {
					fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f v %.1f/%.1f\n",
						here.k, here.x, D.X[D.NX[k] - 1][k], here.v, D.V[D.NV[k] - 1][k]);
					exit(1);
				}

				node_t next = { 0 };
				double u, phi;
				double *Jhere, *Jnext;
				ctrl_t *Chere;
				
				if (here.k == 0)
					assert(here.ix == 0 && here.iv == 0);

				/* final stage costs */
				if (here.k == P.numsteps - 1) {
					J[here.k][here.ix][here.iv] = J_opt[here.k][here.ix][here.iv];
					continue;
				}
				
				/* intermediate state */
				Jhere = &J[here.k][here.ix][here.iv];
				Chere = &C[here.k][here.ix][here.iv];
			
				*Jhere = DBL_MAX;
				next.k = here.k + 1;

				/* for each possible u */
				for (u = D.umin; u <= D.umax; u += D.du) {
					next.x = here.x + here.v * P.T + 0.5 * u * pow(P.T, 2);
					next.v = here.v + u * P.T;
						
					if (discretizexv(&next, next.k) != 0)
						continue;			

					if (crash(next) > 0 && OBSTACLES)
						continue;

					Jnext = &J[next.k][next.ix][next.iv];
					phi = ((1.0 - p[here.k]) * (0.5*pow(u, 2) + *Jnext)) + (p[here.k] * J_opt[here.k][here.ix][here.iv]);
					
					if (phi < *Jhere) {
						*Jhere = phi;
						Chere->u = u;
					}
				}
			}
		}
	}
}

static void read_sdp_init() {
	double dval, dval2, dval3;
	int ival;

	FILE *f;
	char buff[0x100];
	snprintf(buff, sizeof(buff), "deterministic/deterministic x0=%.1f - v0=%.1f.txt", P.x0, P.v0);
	f = fopen(buff, "r");

	while (fgets(buff, sizeof(buff), f)) {
		// fprintf(stderr, "%s \n", buff);
		if (buff[0] == '\n')
			break;

		if (sscanf(buff, "%d %lf %lf %lf", &ival, &dval, &dval2, &dval3) == 4) {
			// fprintf(stderr, "%.4f \n", dval2);
			initial_x[ival] = dval;
			initial_v[ival] = dval2;
			initial_u[ival] = dval3;
		}

	}

}

int main(int argc, char **argv) {
	int k;
	
	P.numsteps = (int)((T_max - T_tl)) + 1;
	P.T = TIME_STEP;
  	// du = init_du;

  	double scenario_x[] = { 0.0,  10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 };
  	double scenario_v[] = { 1.0,   2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0, 11.0 };

	// double loop for different initial conditions
	for (int xi = 0; xi < 1; xi++) {
		for (int vi = 0; vi < 1; vi++) {

      du = init_du;
      // P.x0 = scenario_x[xi]; P.v0 = scenario_v[vi];
      P.x0 = 50.0; P.v0 = 11.0;
			fprintf(stderr, "x0: %.4f -- v0: %.4f \n", P.x0, P.v0);
      
			assert((opt_x = (double*)calloc(sizeof(double), P.numsteps)));
			assert((opt_v = (double*)calloc(sizeof(double), P.numsteps)));
			assert((opt_u = (double*)calloc(sizeof(double), P.numsteps)));

			assert((init_x = (double*)calloc(sizeof(double), P.numsteps)));
			assert((init_v = (double*)calloc(sizeof(double), P.numsteps)));
			assert((init_u = (double*)calloc(sizeof(double), P.numsteps)));

			assert((initial_x = (double*)calloc(sizeof(double), P.numsteps)));
			assert((initial_v = (double*)calloc(sizeof(double), P.numsteps)));
			assert((initial_u = (double*)calloc(sizeof(double), P.numsteps)));
      
			read_sdp_init();

			// for (int k = 0; k < P.numsteps; k++)
			// 	fprintf(stderr, "(%d) \t init_x: %.4f \t init_v: %.4f \t init_a: %.4f \n", k, initial_x[k], initial_v[k], initial_u[k]);

			for (k = 0; k < P.numsteps; k++) {
				opt_x[k] = initial_x[k];
				opt_v[k] = initial_v[k];
				opt_u[k] = initial_u[k];
                // opt_x[k] = in_x[k];
				// opt_v[k] = in_v[k];
				// opt_u[k] = in_a[k];
			}

			node_t initial = { 0 };
			initial.k = 0;
			initial.x = P.x0;
			initial.v = P.v0;
			initial.ix = initial.iv = 0;
			
			clock_t start_all, end_all;
			start_all = clock();
			double opt_J_prev_it = DBL_MAX, du_it_prev = du;
			int it = 0;
			while(true){
				
				for (k = 0; k < P.numsteps; k++) {
                    init_x[k] = opt_x[k];
                    init_v[k] = opt_v[k];
                    init_u[k] = opt_u[k];
				}
				
				init_D(P.T, it);
				init_JC();

				/* run dp */
				clock_t start_it, end_it;
				start_it = clock();

				dp();

				end_it = clock();
				double cpu_time_used_it = ((double)(end_it - start_it)) / CLOCKS_PER_SEC;
				
				/*  prints */
				printsol(initial, it);
				
				if (DISPLAY) {
                    fprintf(stderr, "\n--------------------------------\n");
                    fprintf(stderr, "DU: %.2f -- DV: %.2f DX: %.2f \n", D.du, D.dv, D.dx);
                    fprintf(stderr, "Iteration %d CPU time: %.4f \n", it + 1, cpu_time_used_it);
                    fprintf(stderr, "--------------------------------\n");
				}
				else
				    fprintf(stderr, "Iteration %d -- DU=DV: %.5f -- DX: %.5f -- cost: %.5f -- CPU time: %.4f \n", it + 1, D.du, D.dx, opt_J_it, cpu_time_used_it);
			
				/* terminal criterion: || u(k) - u(k-1) || < epsilon */
				double sum = 0.0;
				for (int i = 0; i < P.numsteps; i++)
				    sum += pow(opt_u[i] - init_u[i], 2);
				double term_criterion = sqrt(sum);

				/* check if reduction of du is needed */
				du_it_prev = du;
				if (opt_J_it == opt_J_prev_it)
				du /= 2.0;
				
				free_JC();
				
				if (opt_J_it > opt_J_prev_it) {
                    fprintf(stderr, "Termination error: current cost > previous cost !!! \n");
                    break;
				}
				if (TERM_CRITERION_DU && (du < 0.125)) {
                    fprintf(stderr, "Termination criterion: DU < 0.125 \n");
                    break;
				}

				opt_J_prev_it = opt_J_it;
				it++;
				
			}
			end_all = clock();
			double cpu_time_used_all = ((double)(end_all - start_all)) / CLOCKS_PER_SEC;
			fprintf(stderr, "Total CPU time (in %d it.): %.4f \n", it + 1, cpu_time_used_all);


			free(opt_x); free(opt_v); free(opt_u);
			free(init_x); free(init_v); free(init_u);
			free(initial_x); free(initial_v); free(initial_u);

			for (int i = 0; i < P.numsteps; i++) {
				free(D.V[i]);
				free(D.X[i]);
			}
			free(D.V);
			free(D.X);

			free(D.vmax); free(D.vmin); free(D.xmax); free(D.xmin);
			free(D.NV); free(D.NX);
			free(D.NXlive_min); free(D.NXlive_max); free(D.NVlive_min); free(D.NVlive_max);
			
			fprintf(stderr, " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n ");
		}
	}


	return 0;
}
