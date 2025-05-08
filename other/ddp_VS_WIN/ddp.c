#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdbool.h>

/* scenario 1 */
#include "../../results_ddp_mathematica/coefficients/coefficients_x0_0_v0_11.h"
//#include "coefficients.h"
double X0 = 0.0, V0 = 11.0;
double const_off = 0.1;

#define Power(a, b) pow(a, b)
//#define sscanf_s sscanf
#define MAX(a,b) ((a > b)?(a):(b))
#define MIN(a,b) (a < b)?(a):(b)
#define ACCESS(M, state) M[state.k][state.ix][state.iv]

#define STOCHASTIC 1
#define TRIANGULAR 0
#define FIXED_TIME 0
#define TIME_STEP 1.0
#define TERM_CRITERION 0.1
#define TERM_CRITERION_2 0.00001



double T_min = 10.0, T_max = 30.0, T_tl = 0.0;  	// possible change times from T_min to T_max

/* analytic parameters  */
double ve = 11.0, xe = 220.0;						// analytic optimal control final states
double w = 0.1;
double tf = 38.0;
double wopt = 1.0;

double *init_x, *init_v, *init_u;
double *x_opt, *v_opt, *u_opt;
/* UP problem */
double te_up;

#define OUTPUT 1
//#define SCENARIO 2

double *initial_v;
double *initial_x;
double *initial_a;

double *phi;

/* problem */
struct {
	double T;
	double x0, y0, v0;
	int numsteps;
}P = { 0 };

struct ddp {
	double **D, **E, **F, **Finv, *G, *H;
	double *alpha, **beta;
	double **A, *B, c;
	double **temp, **temp2;
	double *tempOneDim, *tempOneDim2;
}stage[31];

static void free_JC(void) {
	int k, i;
	for (k = 0; k < P.numsteps; k++) {
		for (i = 0; i < 2; i++) {
			free(stage[k].temp[i]);
			free(stage[k].temp2[i]);
			free(stage[k].D[i]);
			free(stage[k].E[i]);
			free(stage[k].F[i]);
			free(stage[k].Finv[i]);

			free(stage[k].beta[i]);
			free(stage[k].A[i]);
		}
		free(stage[k].temp);
		free(stage[k].temp2);
		free(stage[k].D);
		free(stage[k].E);
		free(stage[k].F);
		free(stage[k].Finv);
		free(stage[k].A);
		free(stage[k].beta);

		free(stage[k].G);
		free(stage[k].H);
		free(stage[k].alpha);
		free(stage[k].B);
		free(stage[k].tempOneDim);
		free(stage[k].tempOneDim2);
	}
	/*free(x_opt);
	free(v_opt);
	free(u_opt);*/

}

static void mat_mult(double **mat1, double **mat2, double **res, int dim1)
{
    int i, j, k;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
            res[i][j] = 0;
            for (k = 0; k < dim1; k++)
                res[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
}

static void mat_mult_vect(double **mat1, double *mat2, double *res, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
            res[i] += mat1[i][j] * mat2[j];
        }
    }
}

static void vect_mult_mat(double *mat1, double **mat2, double *res, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
            res[i] += mat1[j] * mat2[j][i];
        }
    }
}  

static void mat_mult_const(double **mat1, double mat2, double **res, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
                res[i][j] = mat1[i][j] * mat2;
        }
    }
}

static void vect_mult_const(double *mat1, double mat2, double *res, int dim1)
{
    int i;
    for (i = 0; i < dim1; i++) {
		res[i] = mat1[i] * mat2;
    }
}

static void mat_sub(double **mat1, double **mat2, double **res, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
                res[i][j] = mat1[i][j] - mat2[i][j];
        }
    }
}

static void mat_sum(double **mat1, double **mat2, double **res, int dim1)
{
    int i, j;
    for (i = 0; i < dim1; i++) {
        for (j = 0; j < dim1; j++) {
                res[i][j] = mat1[i][j] + mat2[i][j];
        }
    }
}

static void mat_sub2(double *mat1, double *mat2, double *res, int dim1)
{
    int i;
    for (i = 0; i < dim1; i++) {
		res[i] = mat1[i] - mat2[i];
    }
}

static void mat_sum2(double *mat1, double *mat2, double *res, int dim1)
{
    int i;
    for (i = 0; i < dim1; i++) {
		res[i] = mat1[i] + mat2[i];
    }
}

static void mat_2x2_inv(double **mat1, double **res)
{
    
	double det = (mat1[0][0] * mat1[1][1]) - (mat1[0][1] * mat1[1][0]);
	if (det == 0)
		fprintf(stderr, "Determinant equal to zero \n");
	
    res[0][0] = (1.0/det) * mat1[1][1];
	res[0][1] = (1.0/det) * (-mat1[0][1]);
	res[1][0] = (1.0/det) * (-mat1[1][0]);
	res[1][1] = (1.0/det) * mat1[0][0];
}

static void mat_2x2_transpose(double **mat, double **res)
{
    
	int i, j;
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++)
    		res[j][i] = mat[i][j];
	}
}

static double vect_mult_vect(double *mat1, double *mat2, int dim1)
{
	double res = 0.0;

	for (int i = 0; i < dim1; i++) {
		res += mat1[i] * mat2[i];
	}

	return res;
}

static void vect_mult_vect_to_mat(double *mat1, double *mat2, double **res, int dim1)
{

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim1; j++) {
			res[i][j] = mat1[i] * mat2[j];
		}
	}

}

double constraint(double dx, double dv, double xp, double vp, double ap) {
	double res;
	double xU = 150.0;

	return res = -(2.0 * (xU - xp - vp * P.T - 0.5 * ap * pow(P.T, 2)) / pow(P.T, 2)) + 2.0 * dx/pow(P.T, 2) + 2.0 * dv/pow(P.T, 2);
}

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
static double UP_control(double t, double t0, double x0, double v0, double te) {
		
	double res = (5.999999999999999*t*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)));

	return res;
}

static double te_func(double t0, double x0, double v0, double te)
{
	double res = 0.5*w + 0.5*Power((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 +
		2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
			5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
			3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
			3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
			3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))), 2) +
		((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
				5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
				3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
				3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
				3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((-5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 +
			2.*t0*xe - 2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
				5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
				3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
				3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
				3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
		(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)*((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve -
				2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
				(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 +
					5.000000000000001*t0*Power(te, 3)*v0 - 2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve +
					3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve + 0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve -
					3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 + 3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 +
					3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe - 3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
				(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 +
					1.0000000000000004*t0*Power(te, 4)*v0 - 1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve -
					2.999999999999999*Power(t0, 3)*Power(te, 2)*ve + 4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve -
					5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 - 5.999999999999999*t0*Power(te, 3)*x0 +
					5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4));

	return res;
}

static double dte_func(double t0, double x0, double v0, double te)
{
	
	double res = 1.*((5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
		(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
			5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
			3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
			9.000000000000002*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
				2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
				0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
				3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
				3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
			((-5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
				(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
				(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) -
				(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
				(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
					5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
					3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
					9.000000000000002*Power(te, 2)*xe)) /
					((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
				(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
				(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
					2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
					0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
					3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
					3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) +
				(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
					2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
					0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
					3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
					3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
					(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
				((5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
					(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
						2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
						0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
						3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
						3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
						((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
					((5.999999999999999*te*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
						(5.999999999999999*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
						(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
						Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
						(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
						(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
						(1.9999999999999996*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
							5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
							3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
							9.000000000000002*Power(te, 2)*xe)) /
							((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
		((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
		(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))*
		((-5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)))) +
			(5.999999999999999*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)*
		((2.9999999999999996*Power(te, 2)*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)) /
				(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(2.9999999999999996*Power(te, 2)*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2) +
		(5.999999999999999*te*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*te*(-0.9999999999999998*Power(t0, 3)*v0 - 6.000000000000001*Power(t0, 2)*te*v0 + 15.000000000000004*t0*Power(te, 2)*v0 - 8.000000000000002*Power(te, 3)*v0 -
			5.000000000000001*Power(t0, 3)*ve + 6.000000000000001*Power(t0, 2)*te*ve + 2.999999999999999*t0*Power(te, 2)*ve - 3.9999999999999996*Power(te, 3)*ve +
			3.0000000000000004*Power(t0, 2)*x0 + 6.000000000000001*t0*te*x0 - 9.000000000000002*Power(te, 2)*x0 - 3.0000000000000004*Power(t0, 2)*xe - 6.000000000000001*t0*te*xe +
			9.000000000000002*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.*(1.9999999999999996*Power(t0, 4)*v0 - 9.999999999999996*Power(t0, 3)*te*v0 + 8.999999999999996*Power(t0, 2)*Power(te, 2)*v0 + 4.000000000000002*t0*Power(te, 3)*v0 -
			5.000000000000001*Power(te, 4)*v0 - 1.0000000000000004*Power(t0, 4)*ve - 5.999999999999998*Power(t0, 3)*te*ve + 14.999999999999995*Power(t0, 2)*Power(te, 2)*ve -
			7.999999999999997*t0*Power(te, 3)*ve - 5.999999999999999*Power(t0, 3)*x0 + 23.999999999999996*Power(t0, 2)*te*x0 - 17.999999999999996*t0*Power(te, 2)*x0 +
			5.999999999999999*Power(t0, 3)*xe - 23.999999999999996*Power(t0, 2)*te*xe + 17.999999999999996*t0*Power(te, 2)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.9999999999999996*te*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) -
		(1.9999999999999996*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
		(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) -
		(1.*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2)) +
		(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			(Power(1.*t0 - 1.*te, 2)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) +
		(5.999999999999999*(-2.*t0*v0 + 2.*te*v0 - 2.*t0*ve + 2.*te*ve + 2.*x0 - 2.*xe)*
		((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
		2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
			2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
			0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
			3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
			3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
		(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
			1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
			4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
			5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
			((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
		(5.999999999999999*(-4.*Power(t0, 3) + 11.999999999999998*Power(t0, 2)*te - 12.*t0*Power(te, 2) + 4.000000000000001*Power(te, 3))*
		(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe - 2.*te*xe)*
		((2.9999999999999996*Power(te, 2)*(1.*Power(t0, 2)*v0 - 2.*t0*te*v0 + 1.*Power(te, 2)*v0 + 1.*Power(t0, 2)*ve - 2.*t0*te*ve + 1.*Power(te, 2)*ve - 2.*t0*x0 + 2.*te*x0 + 2.*t0*xe -
			2.*te*xe)) / (1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4)) -
			(1.9999999999999996*te*(1.*Power(t0, 4)*v0 - 0.9999999999999998*Power(t0, 3)*te*v0 - 3.0000000000000004*Power(t0, 2)*Power(te, 2)*v0 + 5.000000000000001*t0*Power(te, 3)*v0 -
				2.0000000000000004*Power(te, 4)*v0 + 2.0000000000000004*Power(t0, 4)*ve - 5.000000000000001*Power(t0, 3)*te*ve + 3.0000000000000004*Power(t0, 2)*Power(te, 2)*ve +
				0.9999999999999998*t0*Power(te, 3)*ve - 0.9999999999999999*Power(te, 4)*ve - 3.0000000000000004*Power(t0, 3)*x0 + 3.0000000000000004*Power(t0, 2)*te*x0 +
				3.0000000000000004*t0*Power(te, 2)*x0 - 3.0000000000000004*Power(te, 3)*x0 + 3.0000000000000004*Power(t0, 3)*xe - 3.0000000000000004*Power(t0, 2)*te*xe -
				3.0000000000000004*t0*Power(te, 2)*xe + 3.0000000000000004*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))) +
			(1.*(1.9999999999999996*Power(t0, 4)*te*v0 - 4.999999999999998*Power(t0, 3)*Power(te, 2)*v0 + 2.999999999999999*Power(t0, 2)*Power(te, 3)*v0 + 1.0000000000000004*t0*Power(te, 4)*v0 -
				1.0000000000000002*Power(te, 5)*v0 + 1.*Power(t0, 5)*ve - 1.0000000000000004*Power(t0, 4)*te*ve - 2.999999999999999*Power(t0, 3)*Power(te, 2)*ve +
				4.999999999999998*Power(t0, 2)*Power(te, 3)*ve - 1.9999999999999993*t0*Power(te, 4)*ve - 5.999999999999999*Power(t0, 3)*te*x0 + 11.999999999999998*Power(t0, 2)*Power(te, 2)*x0 -
				5.999999999999999*t0*Power(te, 3)*x0 + 5.999999999999999*Power(t0, 3)*te*xe - 11.999999999999998*Power(t0, 2)*Power(te, 2)*xe + 5.999999999999999*t0*Power(te, 3)*xe)) /
				((1.*t0 - 1.*te)*(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4))))) /
		Power(1.*Power(t0, 4) - 4.*Power(t0, 3)*te + 5.999999999999999*Power(t0, 2)*Power(te, 2) - 4.*t0*Power(te, 3) + 1.0000000000000002*Power(te, 4), 2);

			return res;
}

static double newton_raphson(double t0, double x0, double v0, double te) {
	
	/* Newton - Raphson Method */
	int iter = 0;
	double h = te_func(t0, x0, v0, te) / dte_func(t0, x0, v0, te);
	while (fabs(h) >= 0.001 && iter < 100)
	{
		h = te_func(t0, x0, v0, te) / dte_func(t0, x0, v0, te);
		  
		te = te - h;
		iter += 1;
	}

	return te;
}

double compute_real_cost(){

	double J_opt, p;
	double T_max = 30.0;
	int newtonIntervals = 1;
	double newton_add_term = 5.0;
	double test_x[] = {0.00000, 4.93750, 9.75000, 14.43750, 19.00000, 23.43750, 27.81250, 32.18750, 36.56250, 40.93750, 45.31250, 49.68750, 54.06250, 58.43750, 62.81250, 67.18750, 71.62500, 76.18750, 80.87500, 85.68750, 90.62500, 95.68750, 100.87500, 106.18750, 111.68750, 117.43750, 123.43750, 129.68750, 136.18750, 142.93750, 150.00000};
	double test_v[] = {5.00000, 4.87500, 4.75000, 4.62500, 4.50000, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.37500, 4.50000, 4.62500, 4.75000, 4.87500, 5.00000, 5.12500, 5.25000, 5.37500, 5.62500, 5.87500, 6.12500, 6.37500, 6.62500, 6.87500, 7.25000};
	double test_a[] = {-0.12500, -0.12500, -0.12500, -0.12500, -0.12500, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.12500, 0.12500, 0.12500, 0.12500, 0.12500, 0.12500, 0.12500, 0.12500, 0.25000, 0.25000, 0.25000, 0.25000, 0.25000, 0.25000, 0.37500, 0.00000};

	phi = (double*)calloc(P.numsteps, sizeof(double));
	for (int k = P.numsteps-1; k >= 0; k--) {
		
		double x = x_opt[k], v = v_opt[k], u = u_opt[k];

		double result = 0.0;
		if (T_min <= k && k <= T_max) {
			double te_init, temp_te, check_te, temp_cost = DBL_MAX;
			double checkSpeed = v;
			int negativeSpeed = 0;
			
			
			/* calculate max te assuming the min speed is 1 */
			double max_te = (xe - x) / 1.0;
			te_init = k + 1.0;

			for (int i = 0; i < newtonIntervals; i++) {
				double te_cost = 0.0;

				check_te = newton_raphson(k, x, v, te_init);

				if (check_te > 2.0*max_te || isnan(check_te) || check_te < 0.0)
					continue;

				for (int i = 0; i < P.numsteps; i++) {
					double step = (check_te - k) / (P.numsteps);
					double t = (i*(check_te - k) / (P.numsteps - 1.0) + k);

					te_cost += 0.5*pow(UP_control(t, k, x, v, check_te), 2);
					/* skip optimal te that leads to negative speed trajectories */
					checkSpeed += UP_control(t, k, x, v, check_te)*P.T;
					if (checkSpeed < 0.0)
						negativeSpeed = 1;
				}

				if ((check_te > k) && (te_cost < temp_cost) && (negativeSpeed != 1)) {
					temp_cost = te_cost;
					temp_te = check_te;
				}
				te_init += newton_add_term;
			}
			double te = temp_te;			
			double te_opt = te;
			
			// calculate the optimal control cost
			double step = 0.0;
			for (int i = 0; i < P.numsteps; i++) {
				step = (te - k) / (P.numsteps);
				double t = (i*(te - k) / (P.numsteps - 1.0) + k);
				result += 0.5* (pow(UP_control(t, k, x, v, te), 2));
			}
			result *= step;
			
			J_opt = wopt * result;
			p = 1.0 / ((T_max - k) + 1.0);
		}
		else 
			p = 0.0;

		if (k == P.numsteps - 1) {
			phi[k] = J_opt;
			continue;
		}
		// fprintf(stderr, "(%d) phi: %.4f -- p: %.4f -- J_opt: %.4f \n", k, phi[k+1], p, J_opt);
		
		phi[k] = ((1.0 - p) * (0.5*pow(u, 2) + phi[k+1])) + (p * J_opt);
		fprintf(stderr, "(%d) phi: %.4f \n", k, phi[k]);
	}
	return 0;
}

static void dp(int it) {
	
	int k;
	
	x_opt = (double*)calloc(P.numsteps, sizeof(double));
	v_opt = (double*)calloc(P.numsteps, sizeof(double));
	u_opt = (double*)calloc(P.numsteps, sizeof(double));

	for (k = P.numsteps - 1; k >= 0; k--) {
		((stage[k].temp = (double **)calloc(sizeof(double*), 2)));
		((stage[k].temp2 = (double **)calloc(sizeof(double*), 2)));
		((stage[k].D = (double **)calloc(sizeof(double*), 2)));
		((stage[k].E = (double **)calloc(sizeof(double*), 2)));
		((stage[k].F = (double **)calloc(sizeof(double*), 2)));
		((stage[k].Finv = (double **)calloc(sizeof(double*), 2)));
		((stage[k].A = (double **)calloc(sizeof(double*), 2)));
		((stage[k].beta = (double **)calloc(sizeof(double*), 2)));
		for (int i = 0; i < 2; i++) {
			((stage[k].temp[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].temp2[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].D[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].E[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].F[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].Finv[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].beta[i] = (double *)calloc(sizeof(double), 2)));
			((stage[k].A[i] = (double *)calloc(sizeof(double), 2)));
		}
		((stage[k].G = (double *)calloc(sizeof(double), 2)));
		((stage[k].H = (double *)calloc(sizeof(double), 2)));
		((stage[k].alpha = (double *)calloc(sizeof(double), 2)));
		((stage[k].B = (double *)calloc(sizeof(double), 2)));
		((stage[k].tempOneDim = (double *)calloc(sizeof(double), 2)));
		((stage[k].tempOneDim2 = (double *)calloc(sizeof(double), 2)));
		
		/* fitting function */
		int fitCount = 0;
		double p1 = 0.0, p2 = 0.0, p3 = 0.0, p4 = 0.0, p5 = 0.0, p6 = 0.0;
		double prob = 0;
	
		if (k >= T_min) {
			fitCount = (P.numsteps - 1) - k;

			p1 = fitCoeff[it][fitCount][0]; p2 = fitCoeff[it][fitCount][1]; p3 = fitCoeff[it][fitCount][2];
			p4 = fitCoeff[it][fitCount][3]; p5 = fitCoeff[it][fitCount][4]; p6 = fitCoeff[it][fitCount][5];
			
			prob = 1.0 / ((P.numsteps - 1.0 - k) + 1.0);
			// fprintf(stderr, "p1: %.6f -- p2: %.6f -- p3: %.6f -- p4: %.6f -- p5: %.6f -- p6: %.6f \n", p1, p2, p3, p4, p5, p6);
		}

		double xPoint = init_x[k], vPoint = init_v[k], aPoint = init_u[k];

		if (k == P.numsteps - 1) {

			// create D[K-1] base on the fitting function
			stage[k].D[0][0] = 2.0 * p1;
			stage[k].D[0][1] = p3;
			stage[k].D[1][0] = p3;
			stage[k].D[1][1] = 2.0 * p2;

			// create G[K-1] base on the fitting function
			stage[k].G[0] = p4;
			stage[k].G[1] = p5;

			// for K-1 all the other vectors and matrices are zero, as u is zero
			// so alpha[K-1] and beta[K-1] are zero vectores

			// V[x(K-1), K-1] = 1/2 A(K-1) dx(K-1)^2 + B(K-1) dx(K-1) + (constant)
			// where A(K-1) = D(K-1) - E(K-1) F(K-1)^-1 E(K-1)
			// and   B(K-1) = G(K-1) - E(K-1) F(K-1)^-1 H(K-1)

			stage[k].A[0][0] = stage[k].D[0][0];
			stage[k].A[0][1] = stage[k].D[0][1];
			stage[k].A[1][0] = stage[k].D[1][0];
			stage[k].A[1][1] = stage[k].D[1][1];

			stage[k].B[0] = stage[k].G[0];
			stage[k].B[1] = stage[k].G[1];

			stage[k].c = p6;

			continue;
		}
		
		// define F
		stage[k].F[0][0] = (1.0 - prob) * (1.0 + 0.5 * (0.0 + 1.0 * pow(P.T, 2) * (stage[k+1].A[1][0] * P.T + 0.5 * stage[k+1].A[0][0] * pow(P.T, 2)) + 2.0 * P.T * (stage[k+1].A[1][1] * P.T + 0.5 * stage[k+1].A[0][1] * pow(P.T, 2))));
		stage[k].F[0][1] = 0.0;
		stage[k].F[1][0] = 0.0;
		stage[k].F[1][1] = 1.0;
		
		// define H
		stage[k].H[0] = aPoint * (1.0 - prob) + stage[k+1].B[1] * (1.0 - prob) * P.T + 0.5 * stage[k+1].B[0] * (1.0 - prob) * pow(P.T, 2);
		stage[k].H[1] = 0.0;

		// define D
		stage[k].D[0][0] = stage[k+1].A[0][0] * (1.0 - prob) + 2.0 * p1 * prob;
		stage[k].D[0][1] = 0.5 * (2 * stage[k+1].A[0][0] + stage[k+1].A[0][1] + stage[k+1].A[1][0]) * (1. - prob) + p3 * prob;
		stage[k].D[1][0] = 0.5 * (2 * stage[k+1].A[0][0] + stage[k+1].A[0][1] + stage[k+1].A[1][0]) * (1. - prob) + p3 * prob;
		stage[k].D[1][1] = 0.5 * (2 * (stage[k+1].A[0][0] + stage[k+1].A[1][0]) + 2.0 * (stage[k+1].A[0][1] + stage[k+1].A[1][1])) * (1. - prob) + 2.0 * p2 * prob;
		
		// define E
		stage[k].E[0][0] = 0.5 * stage[k+1].A[0][1] * (1.0 - prob) * P.T + 0.5 * stage[k+1].A[1][0] * (1.0 - prob) * P.T + 0.5 * stage[k+1].A[0][0] * (1.0 - prob) * pow(P.T, 2);
		stage[k].E[0][1] = 0.0;
		stage[k].E[1][0] = 0.5 * stage[k+1].A[0][1] * (1.0 - prob) * P.T + 0.5 * stage[k+1].A[1][0] * (1.0 - prob) * P.T + 1.0 * stage[k+1].A[1][1] * (1.0 - prob) * P.T + 0.5 * stage[k+1].A[0][0] * (1.0 - prob) * pow(P.T, 2) + 0.25 * stage[k+1].A[0][1] * (1. - prob) * pow(P.T, 2) + 0.25 * stage[k+1].A[1][0] * (1.0 - prob) * pow(P.T, 2);
		stage[k].E[1][1] = 0.0;

		// define G
		stage[k].G[0] = stage[k+1].B[0] * (1.0 - prob) + p4 * prob;
		stage[k].G[1] = stage[k+1].B[0] * (1.0 - prob) + stage[k+1].B[1] * (1.0 - prob) + p5 * prob;

		double xvPoint[2];
		xvPoint[0] = xPoint;xvPoint[1] = vPoint;
		double qq[2];
		// qq = -HH - EE . ( {{xPoint},{vPoint}} );
		mat_mult_vect(stage[k].E, xvPoint, stage[k].tempOneDim, 2);
		vect_mult_const(stage[k].H, -1, stage[k].tempOneDim2, 2);
		mat_sub2(stage[k].tempOneDim2, stage[k].tempOneDim, qq, 2);
		
		double *Ct;
		((Ct = (double *)calloc(sizeof(double), 2)));
		// Ct = ( {{0.0},{0.0}} ); d = {{0.0}};
		Ct[0] = 0.0; Ct[1] = 0.0;
		double d = 0.0;
		double constr = constraint(0.0, 0.0, xPoint, vPoint, aPoint);		
		if (constr + const_off >= -aPoint) {
			Ct[0] = -1.0; Ct[1] = 0.0;
			d = constr;
		}
		
		/* Quadratic Problem for Active Constraints */
		double tempCstar;
		stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;

		// tempCstar = Transpose[Ct] . Inverse[FF] . Ct;
		mat_2x2_inv(stage[k].F, stage[k].Finv);
		mat_mult_vect(stage[k].Finv, Ct, stage[k].tempOneDim, 2);
		tempCstar = vect_mult_vect(Ct, stage[k].tempOneDim, 2);
		
		if (tempCstar != 0.0)
			tempCstar = 1.0 / tempCstar;
		// else 
		// 	fprintf(stderr, "Cstar is zero \n");

		double *Cstar, **Hstar;
		Cstar = (double *)calloc(2, sizeof(double));
		Hstar = (double **)calloc(2, sizeof(double*));
		for (int i = 0; i < 2; i++)
			Hstar[i] = (double *)calloc(2, sizeof(double));
		
		stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;

		// Cstar = tempCstar . Transpose[Ct] . Inverse[FF];
		vect_mult_mat(Ct, stage[k].Finv, stage[k].tempOneDim, 2);
		vect_mult_const(stage[k].tempOneDim, tempCstar, Cstar, 2);

		// Hstar = Inverse[FF] . (( {{1.0, 0.0},{0.0, 1.0}} ) - Ct . Cstar);
		double **eyeMat;
		eyeMat = (double **)calloc(2, sizeof(double*));
		for (int i = 0; i < 2; i++)
			eyeMat[i] = (double *)calloc(2, sizeof(double));

		eyeMat[0][0] = 1.0; eyeMat[0][1] = 0.0; eyeMat[1][0] = 0.0; eyeMat[1][1] = 1.0;
		stage[k].temp[0][0] = 0.0; stage[k].temp[0][1] = 0.0; stage[k].temp[1][0] = 0.0; stage[k].temp[1][1] = 0.0;
		vect_mult_vect_to_mat(Ct, Cstar, stage[k].temp, 2);
		mat_sub(eyeMat, stage[k].temp, stage[k].temp2, 2);
		mat_mult(stage[k].Finv, stage[k].temp2, Hstar, 2);
		
		// pStar = Hstar . qq + Transpose[Cstar] . d;

		stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;
		// lambdaStar = Cstar . qq - tempCstar . d;
		double temp = vect_mult_vect(Cstar, qq, 2);
		double lambdaStar = temp - tempCstar * d;
		
		if (lambdaStar > 0) {
			Ct[0] = -1.0; Ct[1] = 0.0;
			d = constraint(0.0, 0.0, xPoint, vPoint, aPoint);

			//ADD Cstar and Hstar
			mat_2x2_inv(stage[k].F, stage[k].Finv);
			mat_mult_vect(stage[k].Finv, Ct, stage[k].tempOneDim, 2);
			tempCstar = vect_mult_vect(Ct, stage[k].tempOneDim, 2);
			
			// If[tempCstar[[1, 1]] != 0, tempCstar = Inverse[tempCstar]];
			if (tempCstar != 0.0)
				tempCstar = 1.0 / tempCstar;
			
			stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;
			// Cstar = tempCstar . Transpose[Ct] . Inverse[FF];
			vect_mult_mat(Ct, stage[k].Finv, stage[k].tempOneDim, 2);
			vect_mult_const(stage[k].tempOneDim, tempCstar, Cstar, 2);

			eyeMat[0][0] = 1.0; eyeMat[0][1] = 0.0; eyeMat[1][0] = 0.0; eyeMat[1][1] = 1.0;
			stage[k].temp[0][0] = 0.0; stage[k].temp[0][1] = 0.0; stage[k].temp[1][0] = 0.0; stage[k].temp[1][1] = 0.0;
			vect_mult_vect_to_mat(Ct, Cstar, stage[k].temp, 2);
			mat_sub(eyeMat, stage[k].temp, stage[k].temp2, 2);
			mat_mult(stage[k].Finv, stage[k].temp2, Hstar, 2);
		}

		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
		double PP = d;
		double R[2];
		R[0] = 2.0; R[1] = 2.0; 

		/////////////////////    DONE UP TO HERE    //////////////////////
		// alpha = -Hstar . HH + Transpose[Cstar] . P;
		stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;
		stage[k].tempOneDim2[0] = 0.0; stage[k].tempOneDim2[1] = 0.0;
		vect_mult_const(Cstar, PP, stage[k].tempOneDim, 2);
		mat_mult_vect(Hstar, stage[k].H, stage[k].tempOneDim2, 2);
		mat_sub2(stage[k].tempOneDim, stage[k].tempOneDim2, stage[k].alpha, 2);

		// beta = -Hstar . Transpose[EE] + Transpose[Cstar] . R;
		double **Etransp;
		Etransp = (double **)calloc(2, sizeof(double*));
		for (int i = 0; i < 2; i++)
			Etransp[i] = (double *)calloc(2, sizeof(double));
		stage[k].temp[0][0] = 0.0; stage[k].temp[0][1] = 0.0; stage[k].temp[1][0] = 0.0; stage[k].temp[1][1] = 0.0;
		stage[k].temp2[0][0] = 0.0; stage[k].temp2[0][1] = 0.0; stage[k].temp2[1][0] = 0.0; stage[k].temp2[1][1] = 0.0;
		vect_mult_vect_to_mat(Cstar, R, stage[k].temp, 2);		

		mat_2x2_transpose(stage[k].E, Etransp);
		mat_mult(Hstar, Etransp, stage[k].temp2, 2);
		mat_sub(stage[k].temp, stage[k].temp2, stage[k].beta, 2);

		/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

		// PP = DD + Transpose[beta] . FF . beta + EE . beta + Transpose[beta] . Transpose[EE]; (PP -> AA)
		double **betaT, **temp3, **temp4, **temp5, **temp6;
		temp3 = (double **)calloc(2, sizeof(double*));
		temp4 = (double**)calloc(2, sizeof(double*));
		temp5 = (double**)calloc(2, sizeof(double*));
		temp6 = (double**)calloc(2, sizeof(double*));
		betaT = (double**)calloc(2, sizeof(double*));
		for (int i = 0; i < 2; i++) {
			temp3[i] = (double *)calloc(2, sizeof(double));
			temp4[i] = (double*)calloc(2, sizeof(double));
			temp5[i] = (double*)calloc(2, sizeof(double));
			temp6[i] = (double*)calloc(2, sizeof(double));
			betaT[i] = (double*)calloc(2, sizeof(double));
		}
		
		stage[k].temp[0][0] = 0.0; stage[k].temp[0][1] = 0.0; stage[k].temp[1][0] = 0.0; stage[k].temp[1][1] = 0.0;
		stage[k].temp2[0][0] = 0.0; stage[k].temp2[0][1] = 0.0; stage[k].temp2[1][0] = 0.0; stage[k].temp2[1][1] = 0.0;
		mat_2x2_transpose(stage[k].beta, betaT);
		mat_mult(betaT, stage[k].F, stage[k].temp, 2);
		mat_mult(stage[k].temp, stage[k].beta, stage[k].temp2, 2); // 2nd term
		mat_mult(stage[k].E, stage[k].beta, stage[k].temp, 2); // 3rd term (den exei midenistei prin to temp!!!!!!)
		mat_mult(betaT, Etransp, temp3, 2); // 4th term
		mat_sum(stage[k].D, stage[k].temp2, temp4, 2);
		mat_sum(temp4, stage[k].temp, temp5, 2);
		mat_sum(temp5, temp3, stage[k].A, 2);

		// QQ = GG + Transpose[beta] . FF . alpha + EE . alpha + Transpose[beta] . HH; (QQ -> BB)
		double *tempOneDim3, *tempOneDim4, *tempOneDim5, *tempOneDim6;
		((tempOneDim3 = (double *)calloc(2, sizeof(double))));
		((tempOneDim4 = (double*)calloc(2, sizeof(double))));
		((tempOneDim5 = (double*)calloc(2, sizeof(double))));
		((tempOneDim6 = (double*)calloc(2, sizeof(double))));
		stage[k].tempOneDim[0] = 0.0; stage[k].tempOneDim[1] = 0.0;
		stage[k].tempOneDim2[0] = 0.0; stage[k].tempOneDim2[1] = 0.0;
		mat_mult(betaT, stage[k].F, stage[k].temp, 2);
		mat_mult_vect(stage[k].temp, stage[k].alpha, stage[k].tempOneDim, 2); // 2nd term
		mat_mult_vect(stage[k].E, stage[k].alpha, stage[k].tempOneDim2, 2); // 3rd term
		mat_mult_vect(betaT, stage[k].H, tempOneDim3, 2); // 4th term

		mat_sum2(stage[k].G, stage[k].tempOneDim, tempOneDim4, 2);
		mat_sum2(tempOneDim4, stage[k].tempOneDim2, tempOneDim5, 2);
		mat_sum2(tempOneDim5, tempOneDim3, stage[k].B, 2);
	}

	double *cost_opt;
	cost_opt = (double*)calloc(P.numsteps, sizeof(double));
	double eps = 1.0;

	x_opt[0] = init_x[0];
	v_opt[0] = init_v[0];

	FILE *f;
	char buf[0x100];
	snprintf(buf, sizeof(buf), "../outputs/(%02d) x0=%.1f - v0=%.1f.txt", it + 1, x_opt[0], v_opt[0]);
	fopen_s(&f, buf,"w");
	// if (OUTPUT) {
	fprintf(f, "k  Position  Speed  Control  Init_x  Init_v  Init_a  SDP_x  SDP_v  SDP_a \n");
	// }

	for (int i = 0; i < P.numsteps; i++) {
		double dx = (x_opt[i] - init_x[i]), dv = (v_opt[i] - init_v[i]);

		u_opt[i] = init_u[i] + eps * (stage[i].alpha[0]) + stage[i].beta[0][0] * dx + stage[i].beta[0][1] * dv;
		
		u_opt[i] = MIN(MAX(u_opt[i], 2.0 * (-x_opt[i] - v_opt[i] * P.T) / pow(P.T, 2)), 2.0 * (150.0 - x_opt[i] - v_opt[i] * P.T) / pow(P.T, 2));
		u_opt[i] = MIN(MAX(u_opt[i], (-v_opt[i] / P.T)), ((16.0 - v_opt[i]) / P.T));
		u_opt[i] = MAX(MIN(u_opt[i], 3.0), -3.0);

		x_opt[i + 1] = x_opt[i] + v_opt[i] * P.T + 0.5 * pow(P.T, 2) * u_opt[i];
		v_opt[i + 1] = v_opt[i] + P.T * u_opt[i];

		fprintf(f, "%d  %.4f  %.4f  %.4f  %.4f  %.4f  %.4f \n", 
			i, x_opt[i], v_opt[i], u_opt[i], init_x[i], init_v[i], init_u[i]);
	}

}

void read_init_traj(double x0, double v0) {

	int ival;
	double dval, dval2, dval3;

	initial_x = (double*)calloc(P.numsteps + 1, sizeof(double));
	initial_v = (double*)calloc(P.numsteps + 1, sizeof(double));
	initial_a = (double*)calloc(P.numsteps + 1, sizeof(double));

	FILE* f;
	char buff[0x100];

	snprintf(buff, sizeof(buff), "C:/Users/kmountakis/Desktop/ddp thesis/results_sdp/deterministic/deterministic x0=%.1f - v0=%.1f.txt",x0, v0);
	fopen_s(&f, buff, "r");

	while (fgets(buff, sizeof(buff), f)) {
		if (buff[0] == '\n')
			break;

		if (sscanf_s(buff, "%d %lf %lf %lf", &ival, &dval, &dval2, &dval3) == 4) {
			initial_x[ival] = dval;
			initial_v[ival] = dval2;
			initial_a[ival] = dval3;
		}
	}


}

int main(int argc, char **argv) {
		
	int k;
	
	P.numsteps = (int)((T_max - T_tl)) + 1;
	P.T = TIME_STEP;
	P.x0 = X0;
	P.v0 = V0;

	init_x = (double*)calloc(P.numsteps, sizeof(double));
	init_v = (double*)calloc(P.numsteps, sizeof(double));
	init_u = (double*)calloc(P.numsteps, sizeof(double));
	
	read_init_traj(X0, V0);

	for (k = 0; k < P.numsteps; k++) {
		init_x[k] = initial_x[k];
		init_v[k] = initial_v[k];
		if (k < P.numsteps - 1)
			init_u[k] = initial_a[k];
	}

	int it = 0;
	clock_t start_all, end_all;
	start_all = clock();
	while(true){	
		
		/* run dp */
		dp(it);
		
		
		double sum = 0.0;
		for (int i = 0; i < P.numsteps; i++)
			sum += sqrt(pow(init_u[i] - u_opt[i], 2));
		fprintf(stderr, "it: %d -- norm: %.4f \n", it+1, sum);

		if ( sum < 0.001 || it >= it_mathem_ddp - 1 )
		// if ( sum < 0.001 )
			break;

		for (int i = 0; i < P.numsteps; i++) {
			init_x[i] = x_opt[i];
			init_v[i] = v_opt[i];
			init_u[i] = u_opt[i];
		}

		//free_JC();
		
		it++;
		fprintf(stderr, "~~~~~~~~~~~~~~~~~~~~~~~ \n");
	}
	end_all = clock();
	double cpu_time_used_all = ((double)(end_all - start_all)) / CLOCKS_PER_SEC;
	fprintf(stderr, "\n Total CPU time (in %d it.): %.8f \n\n", it + 1, cpu_time_used_all);

	compute_real_cost();
	
	FILE *fcpu, *fcost;
	if (fopen_s(&fcpu, "c_cpu_times.txt", "a+"));
	if (fopen_s(&fcost, "c_costs.txt", "a+"));
	fprintf(fcpu, "cpu_time(%.1f,%.1f): %.4f\n", X0, V0, cpu_time_used_all);
	fprintf(fcost, "cost(%.1f,%.1f): %.4f\n", X0, V0, phi[0]);
	fclose(fcpu);
	fclose(fcost);

	return 0;
}
