/* only one lane change modification added */
/* crashes are now infeasible transitions */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iy][state.iv]

//int CHRASHED = 0;
//double c_gap = 4.0;

typedef struct node {
	int k;
	int ix, iy, iv;
	double x, y, v;
	double _x, _v;
} node_t;

typedef struct {
	double u;
	int dy;
} ctrl_t;

double ****J;
ctrl_t ****C;

/* problem */
struct {
	double vd, T, C, safety;
	double x0, y0, v0, a0;
	double *obst_x, *obst_y, *obst_v;
	int obst_x_len, obst_y_len, obst_v_len;
	int numlanes, numsteps, n;
	int *skip, nskip;

	int trim;
}P = { 0 };

/* solution */
struct {
	double **obst_x;
	double **obst_y;
	double **obst_v;
}S = { 0 };

/* state variable domains */
struct {
	int NX, *NXlive_min, *NXlive_max;
	double dx, *X;

	int NY;
	double *Y;

	int NV, *NVlive_min, *NVlive_max;
	double vmin, vmax, dv, *V;

	double umin, umax, du;
}D = { 0 };

unsigned interpolate = 0;
unsigned interpolate_oddity = 0;

/* load problem in structure P */
static void read_P(void) {
	char buf[1024];
	double dval;
	int ival;

	assert((P.obst_x = (double *)calloc(sizeof(double), 1024)));
	assert((P.obst_y = (double *)calloc(sizeof(double), 1024)));
	assert((P.obst_v = (double *)calloc(sizeof(double), 1024)));

	while (fgets(buf, sizeof(buf), stdin)) {
		if (buf[0] == '\n')
			break;

		// fprintf(stderr, "%s", buf);
		

		if (sscanf(buf, "\"vd\":%lf", &dval) == 1)
			P.vd = dval;
		if (sscanf(buf, "\"numlanes\":%d", &ival) == 1)
			P.numlanes = ival;
		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			P.numsteps = ival;
		if (sscanf(buf, "\"safety\":%lf\n", &dval) == 1)
			P.safety = 1.0 - dval;
		if (sscanf(buf, "\"C\":%lf", &dval) == 1)
			P.C = dval;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			P.T = dval;
		if (sscanf(buf, "\"x(0)\":%lf", &dval) == 1)
			P.x0 = dval;
		if (sscanf(buf, "\"y(0)\":%lf", &dval) == 1)
			P.y0 = dval;
		if (sscanf(buf, "\"v(0)\":%lf", &dval) == 1)
			P.v0 = dval;
		if (sscanf(buf, "\"a(0)\":%lf\n", &dval) == 1)
			P.a0 = dval;
		
		if (sscanf(buf, "\"obst_x(%d,0)\":%lf\n", &ival, &dval) == 2) {
			P.obst_x[ival] = dval;
			P.obst_x_len = MAX(P.obst_x_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_y(%d,0)\":%lf\n", &ival, &dval) == 2) {
			P.obst_y[ival] = dval;
			P.obst_y_len = MAX(P.obst_y_len, ival + 1);
		}
		if (sscanf(buf, "\"obst_v(%d,0)\":%lf\n", &ival, &dval) == 2) {
			P.obst_v[ival] = dval;
			P.obst_v_len = MAX(P.obst_v_len, ival + 1);
		}

 		memset(buf, 0, sizeof(buf));
 	} 


	if (!(P.obst_x_len == P.obst_y_len && P.obst_x_len == P.obst_v_len)) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	P.n = P.obst_x_len;
	assert((P.skip = (int*)calloc(sizeof(int), P.n)));
}

/* compute obstacle trajectories in structure S */
static void
init_S(int numsteps, double T, int n, double *obst_x, double *obst_y, double *obst_v)
{
	int i, k;

	/* allocate */
	assert((S.obst_x = (double**)calloc(sizeof(double*), numsteps)));
	assert((S.obst_y = (double**)calloc(sizeof(double*), numsteps)));
	assert((S.obst_v = (double**)calloc(sizeof(double*), numsteps)));
	for (k = 0; k < numsteps; k++) {
		assert((S.obst_x[k] = (double*)calloc(sizeof(double), n)));
		assert((S.obst_y[k] = (double*)calloc(sizeof(double), n)));
		assert((S.obst_v[k] = (double*)calloc(sizeof(double), n)));
	}

	/* compute */
	for (k = 0; k < numsteps; k++) {
		if (k == 0) {
			for (i = 0; i < n; i++) {
				S.obst_x[k][i] = obst_x[i];
				S.obst_y[k][i] = obst_y[i];
				S.obst_v[k][i] = obst_v[i];
			}
			continue;
		}
		for (i = 0; i < n; i++) {
			S.obst_x[k][i] = S.obst_x[k - 1][i] + S.obst_v[k - 1][i] * T;
			S.obst_y[k][i] = S.obst_y[k - 1][i];
			S.obst_v[k][i] = S.obst_v[k - 1][i];
		}
	}
}

static node_t *
getnext(node_t here, ctrl_t c);

static int
discretizexv(node_t *state);

static void
init_D(double T)
{
	int ix, iy, iv;

	/* y domain */
	D.NY = P.numlanes;
	assert((D.Y = (double*)calloc(sizeof(double), D.NY)));
	for (iy = 0; iy < D.NY; iy++)
		D.Y[iy] = iy;

	/* v domain: dv is chosen arbitrarily small (i.e. 0.5) */
	D.umin = -4;
	D.umax = +4;

	D.dv = 2.0*P.T;
	//D.dv = 0.5;
	D.vmin = fmod(P.v0, D.dv);
	D.vmax = ceil(MAX(1.1*P.v0, 1.5*P.vd));
	//D.vmax = P.v0 + D.umax * P.numsteps;
	D.NV = (int)ceil(D.vmax / D.dv)+1;
	assert((D.V = (double*)calloc(sizeof(double), D.NV)));
	for (iv = 0; iv < D.NV; iv++) 
		D.V[iv] = D.vmin + D.dv * iv;

	/* u domain: du = dv/T */
	D.du = D.dv/T;
	//D.du = D.dv;

	/* x domain: dx = 0.5 * D.du * T^2 */
	D.dx = 0.5 * D.du * pow(T, 2);
	//D.dx = D.du / 2;
	double xmax = (D.vmax * P.T * P.numsteps);
	D.NX = (int)ceil(xmax/D.dx)+1;
	assert((D.X = (double*)calloc(sizeof(double), D.NX)));
	for (ix = 0; ix < D.NX; ix++) 
		D.X[ix] = P.x0 + D.dx * ix;

	fprintf(stderr, "*\t T %.2f C %.2f dv %.2f du %.2f dx %.2f interpolate %d\n", P.T, P.C, D.dv, D.du, D.dx, interpolate);

	int k, ixmin;
	D.NXlive_min = (int*)calloc(sizeof(int), P.numsteps);
	D.NXlive_max = (int*)calloc(sizeof(int), P.numsteps);

	D.NVlive_min = (int*)calloc(sizeof(int), P.numsteps);
	D.NVlive_max = (int*)calloc(sizeof(int), P.numsteps);
	
	for (k=0; k < P.numsteps; k++) {
		ixmin = (int)round((P.x0 + P.v0*P.T*k + 0.5*D.umin*pow(P.T*k, 2) -D.X[0])/D.dx);
		if (ixmin < 0)
			D.NXlive_min[k] = D.NXlive_min[k-1];
		else
			D.NXlive_min[k] = ixmin;
		D.NXlive_max[k] = (int)round((P.x0 + D.vmax*P.T*k - D.X[0])/D.dx);

		D.NVlive_min[k] = 0;
		D.NVlive_max[k] = D.NV;
	}
}

static void
init_JC(void)
{
	assert((J = (double ****)calloc(sizeof(double*), P.numsteps)));
	assert((C = (ctrl_t ****)calloc(sizeof(ctrl_t*), P.numsteps)));
	int k;
	for (k = 0; k < P.numsteps; k++) {

		assert((J[k] = (double ***)calloc(sizeof(double*), D.NX)));
		assert((C[k] = (ctrl_t ***)calloc(sizeof(ctrl_t*), D.NX)));
		int ix;
		for (ix = 0; ix < D.NX; ix++) {
			assert((J[k][ix] = (double **)calloc(sizeof(double*), D.NY)));
			assert((C[k][ix] = (ctrl_t **)calloc(sizeof(ctrl_t*), D.NY)));
			int iy;
			for (iy = 0; iy < D.NY; iy++) {
				assert((J[k][ix][iy] = (double *)calloc(sizeof(double), D.NV)));
				assert((C[k][ix][iy] = (ctrl_t *)calloc(sizeof(ctrl_t), D.NV)));
			}
		}
	}
}

static void
free_JC(void)
{
	int k, ix, iy;
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < D.NX; ix++) {
			for (iy = 0; iy < D.NY; iy++) {
				free(J[k][ix][iy]);
				free(C[k][ix][iy]);
			}
		}
	}
	for (k = 0; k < P.numsteps; k++) {
		for (ix = 0; ix < D.NX; ix++) {
			free(J[k][ix]);
			free(C[k][ix]);
		}
	}

	for (k = 0; k < P.numsteps; k++) {
		free(J[k]);
		free(C[k]);
	}
}

static int
discretizexv(node_t *state)
{
	int ix, iv;

	ix = (int)round((state->x - D.X[0]) / D.dx);
	iv = (int)round((state->v - D.V[0]) / D.dv);

	if (!(ix >= 0 && ix < D.NX)) 
		return 1;
	
	if (!(iv >= 0 && iv < D.NV))  
		return 2;

	state->x = D.X[(state->ix = ix)];
	state->v = D.V[(state->iv = iv)];

	return 0;
}

static double
crash(node_t state, int dy)
{
	int i; 
	double dx;
	double sgap;
	double cost = 0; 
	double L = 5.0;

 	for (i = 0; i < P.n; i++) {
		if (S.obst_y[state.k][i] != state.y || P.skip[i])
			continue;

		dx = fabs(S.obst_x[state.k][i] - state.x);

		if (S.obst_x[state.k][i] - state.x >= 0)
			sgap = ((-0.5*(pow(state.v - (P.safety)*S.obst_v[state.k][i], 2)) / D.umin) + L);
		else
			sgap = ((-0.5*(pow(S.obst_v[state.k][i] - (P.safety)*state.v, 2)) / D.umin) + L);
			
		if (dx <= sgap )
			cost += 1.0;
		
	}
	return cost;
}

static node_t *
getnext(node_t here, ctrl_t c)
{
	static node_t next;
	int err;

	next.k = here.k + 1;
	next.x = next._x = here.x + here.v*P.T + 0.5*c.u*pow(P.T, 2);
	next.v = next._v = here.v + c.u*P.T;
	if ((err = discretizexv(&next))) 
		return NULL;

	next.y = here.y + c.dy;
	next.iy = (int)next.y;
	if (next.y < 0 || next.y >= D.NY) {
		return NULL;
	}


	return &next;
}

static void
printsol(node_t initial)
{
	node_t *nextptr;
	int k = 0, i;
	node_t *state;

	assert((state = (node_t*)calloc(sizeof(node_t), P.numsteps)));

	puts("var Data = {");

	/* trace sequence of states */
	state[0] = initial;
	ctrl_t c;
	for (k = 0; k < P.numsteps-1; k++) {
		c = ACCESS(C, state[k]);
		fprintf(stderr, "*\t k %d  x %.1f y %.1f v %.1f \t J %.1f \t u %.1f dy %d crash %.1f \n",
			k, state[k].x,
			state[k].y,
			state[k].v,
			ACCESS(J, state[k]),
			c.u, c.dy,
			crash(state[k],0));

		if ((nextptr = getnext(state[k], c)) == NULL) {
			fprintf(stderr, "invalid \n");
			exit(1);
		}
		state[k + 1] = *nextptr;
	}

	/* ego trajectory from squence of states/controls */
	double x = state[0].x;
	double v = state[0].v;
	double u = 0, dy = 0;

	for (k = 0; k < P.numsteps; k++) {
		u = ACCESS(C, state[k]).u;
		dy = ACCESS(C, state[k]).dy;

		for (int in = 0; in < P.trim; in++) {
			printf("\"x(%d)\":%.3f, /*%.2f*/\n", P.trim*k + in,  x, state[k].x);
			printf("\"y(%d)\":%.3f,\n", P.trim*k + in,  state[k].y);
			printf("\"vx(%d)\":%.3f, /*%.1f*/\n", P.trim*k + in,  v, state[k].v);
			printf("\"vy(%d)\":%.3f,\n", P.trim*k + in,  0.0);
			printf("\"ux(%d)\":%.3f,\n", P.trim*k + in, u);
			printf("\"uy(%d)\":%.3f,\n", P.trim*k + in, dy);
		}
		puts("");
		x = x + (v * P.T) + (0.5 * u * pow(P.T, 2));
		v = v + (u * P.T);
	}

	/* obstacle trajectories */
	for (k = 0; k < P.numsteps; k++) {
		for (i = 0; i < P.n; i++) {
			for (int in = 0; in < P.trim; in++) {
				printf("\"obst_x(%d,%d)\":%.3f,\n", i, P.trim*k + in, S.obst_x[k][i]);
				printf("\"obst_y(%d,%d)\":%.3f,\n", i, P.trim*k + in, S.obst_y[k][i]);
				printf("\"obst_v(%d,%d)\":%.3f,\n", i, P.trim*k + in, S.obst_v[k][i]);
				puts("");
			}
		}
		puts("");
	}

	printf("\"n\":%d,\n", P.n - P.nskip);
	printf("\"k\":%d,\n", P.trim*P.numsteps);
	printf("\"numlanes\":%d,\n", P.numlanes);
	printf("\"vd\":%f,\n", P.vd);
	printf("\"Step\":%f,\n", P.T/P.trim);
	printf("\"C\":%f, \n", P.C);
	printf("\"safety\":%f, \n", P.safety);
	//printf("\"c_gap\":%f, \n", c_gap);
	printf("\"a(0)\":%f, \n", P.a0);
	printf("};\n");

	free(state);
}

static void
dp(void)
{
	int k, ix, iy, iv;
	double x, y, v;
	for (k = P.numsteps - 1; k >= 0; k--) {
		for (ix = D.NXlive_min[k]; ix < ((k == 0) ? 1 : D.NXlive_max[k]); ix++) {
			x = ((k == 0) ? P.x0 : D.X[ix]);
			for (iy = 0; iy < ((k == 0) ? 1 : D.NY); iy++) {
				y = ((k == 0) ? P.y0 : D.Y[iy]);
				for (iv = D.NVlive_min[k]; iv < ((k == 0) ? 1 : D.NVlive_max[k]); iv++) {
					v = ((k == 0) ? P.v0 : D.V[iv]);

					node_t here = { 0 };
					here.k = k;
					here.x = x;
					here.y = y;
					here.v = v;
					here.ix = ((k == 0) ? 0 : ix);
					here.iy = ((k == 0) ? 0 : (int)y);
					here.iv = ((k == 0) ? 0 : iv);

					if (here.ix > D.NX || here.iy > D.NY || here.iv > D.NV) {
						fprintf(stderr, "invalid state on the grid: k %d x %.1f/%.1f y %.1f/%.1f v %.1f/%.1f\n",
							here.k, here.x, D.X[D.NX - 1], here.y,
							(double)P.numlanes - 1, here.v, D.V[D.NV - 1]);
						exit(1);
					}

					node_t next = { 0 };
					double dy, u, phi;
					double *Jhere, *Jnext;
					ctrl_t *Chere;

					if (here.k == 0)
						assert(here.ix == 0 && here.iy == 0 && here.iv == 0);

					if (here.k == P.numsteps - 1) {
						/* final stage costs */
						J[here.k][here.ix][here.iy][here.iv] = pow(here.v - P.vd, 2) + crash(here, 0);
						continue;
					}

					/* intermediate state */
					Jhere = &J[here.k][here.ix][here.iy][here.iv];
					Chere = &C[here.k][here.ix][here.iy][here.iv];

					*Jhere = DBL_MAX;
					next.k = here.k + 1;
					for (u = D.umin; u <= D.umax; u += D.du) {
					//for (u = -3.3; u <= 3.3; u += 3.3) {
						next.x = here.x + here.v * P.T + 0.5 * u * pow(P.T, 2);
						next.v = here.v + u * P.T;

						if (discretizexv(&next))
							continue;

						for (dy = -1; dy <= 1; dy++) {
							next.y = here.y + dy;

							/* road boundaries check */
							if (!(next.y >= 0 && next.y < D.NY))
								continue;

							/* infeasible areas check */
							//int infeasible_area = 0;
							//for (int j = 0; j < P.n; j++) {
							//	
							//	/*if (next.x <= S.obst_x[next.k][j])
							//		c_gap = 0.1 * next.v * P.C;
							//	else
							//		c_gap = 0.1 * S.obst_v[next.k][j] * P.C;*/
							//	if ( (fabs(next.x - S.obst_x[next.k][j]) <= c_gap * P.C) && (S.obst_y[next.k][j] == next.y) ) {
							//		infeasible_area = 1;
							//		break;
							//	}
							//}
							//if (infeasible_area == 1)
							//	continue;
							
							if (crash(next, dy) != 0)
								continue;

							/* only one lane change allowed */
							if (here.y != P.y0 && dy != 0)
								continue;

							next.iy = (int)next.y;
							assert(next.y == (double)next.iy);

							Jnext = &J[next.k][next.ix][next.iy][next.iv];

							//phi = pow( (P.vd - here.v)/P.vd, 2) + 5.0*pow(u/3.3, 2) + pow(dy, 2) + here.y;
							phi = pow((P.vd - here.v), 2) + 15.0*pow(u, 2) + pow(dy, 2);// +here.y;

							if (phi + *Jnext < *Jhere) {
								*Jhere = phi + *Jnext;
								Chere->u = u;
								Chere->dy = (int)dy;
							}
						}
					}
				}
			}
		}
	}
}


int main(int argc, char **argv) {
	
	read_P();

	P.trim = (int)(1.0 / P.T);
	P.T *= P.trim;
	P.numsteps /= P.trim;

	clock_t tic, toc;
	double cpu_time;

	tic = clock();

	init_S(P.numsteps, P.T, P.n, P.obst_x, P.obst_y, P.obst_v); 
	init_D(P.T);
	init_JC();

	node_t initial = { 0 };
	initial.k = 0;
	initial.x = P.x0;
	initial.y = P.y0;
	initial.v = P.v0;
	initial.ix = initial.iy = initial.iv = 0;

	dp();

	toc = clock();
	// printsol(initial);
	free_JC();
	
	cpu_time = ((double)(toc - tic)) / CLOCKS_PER_SEC;
	fprintf(stderr, "CPU Time: %.4f \n \n", cpu_time);

	return 0;
}

