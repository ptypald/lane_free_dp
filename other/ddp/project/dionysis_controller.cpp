#include "dionysis_controller.h"


double Dionysis_Feedback::distance(double xi, double xj, double yi, double yj, double var) {
	 double result;
    //result = sqrt((pow((xi - xj),2) + var*pow((yi - yj),2)));
    result = sqrt((xi - xj) * (xi - xj) + p * (yi - yj) * (yi - yj));

    return result;
}

double Dionysis_Feedback::sum(double(*var1), int var2) {
	int var4;
    double result = 0;
    for (var4 = 0; var4 <= var2; var4++) {
        result += var1[var4];
    }

    return result;
}

double Dionysis_Feedback::Vfunc(double var1, int var2, double var3, double var4) {
    double result;
    if (var1 > var2) {
        result = 0;
    }
    else {
        result = (double)(var4 * (-(double)3 * pow(((double)var2 - var1), 2) * (var1 - var3) - pow(((double)var2 - var1), 3))) / (pow((var1 - var3), 2));
    }
    
    return result;
}

double Dionysis_Feedback::func1(double var2, double var1) {
    double result = 0;
    if (var2 <= -var1) {
        result = 0;
    }
    else if (var2 < 0) {
        result = 1 / (2 * var1) * pow((var2 + var1), 2);
    }
    else {
        result = 1 / (2 * var1) * (pow(var1, 2) + 2 * var2 * var1);
    }

    return result;
}

double Dionysis_Feedback::Ufunc(double var1, double var2, double var3, int var4) {
    double result = 0;
    if (var4 == 0) {
        if ((-var2 * sqrt(var3 - 1) / sqrt(var3)) < var1 && var1 < (var2 * sqrt(var3 - 1) / sqrt(var3))) {
            result = 0;
        }
        else {
            result = 4 * qu * (((1 / ((var2 * var2) - (var1 * var1))) - (var3 / (var2 * var2))) * ((1 / ((var2 * var2) - (var1 * var1))) - (var3 / (var2 * var2))) * ((1 / ((var2 * var2) - (var1 * var1))) - (var3 / (var2 * var2)))) * ((2 * var1) / (((var2 * var2) - (var1 * var1)) * ((var2 * var2) - (var1 * var1))));
        }
    }
    else
    {
        if (var1 > -var2 + 1 && var1 < var2 - 1)
        {
            result = 0;
        }
        else {
            result = (4 * pow(log((var2 - var1)), 3) * pow(log((var2 + var1)), 3) * ((-var2 + var1) * log((var2 - var1)) + (var2 + var1) * log(var2 + var1))) / (-pow(var2, 2) + pow(var1, 2));
        }
    }
    return result;
}

double Dionysis_Feedback::kappa(double var1, double var2, double var3, double var4) {
    double result = 0;
    if ((var2 <= var1) && (var1 <= var3))
    {
        result = var4 * pow((var3 - var1), 2);
    }
    else
    {
        result = 0;
    }

    return result;
}

int Dionysis_Feedback::mod(int v1, int v2) {
	int ret = v1 % v2;
    if (ret < 0)
        ret += v2;
    
	return ret;
}

void Dionysis_Feedback::run() {
    double *init;
    init = (double*)calloc(4, sizeof(double));

	init[0] = x0;
	init[1] = y0 - a;
	init[2] = 0.0;	// theta
	init[3] = vx0;

    double *t;

    clock_t start, end;
    double cpu_time_used;

    double delta;
    int  i, j, z;
    N = (int)(tmax / T);
    int N2 = N / every + 1;
    double dis[n][n], Vdot[n][n];
    double *V1, *V2, *viscx, *viscy, *viscLx, *viscLmax;

    V1 = (double*)calloc(n, sizeof(double));
    V2 = (double*)calloc(n, sizeof(double));
    viscx = (double*)calloc(n, sizeof(double));
    viscy = (double*)calloc(n, sizeof(double));
    viscLx = (double*)calloc(n, sizeof(double));
    viscLmax = (double*)calloc(n, sizeof(double));


    double u = 0, F = 0, k = 0, Li = 0;
    t = (double*)calloc(N2 + 1, sizeof(double));

	double x = init[0], y = init[1], theta = init[2], v = init[3];
	vector<double> obsx, obsy, obsv, obstheta;
    int obs_n = obsx0.size();

	for (int i = 0; i < obs_n; i++) {
		obsx.push_back(obsx0[i]);
		obsy.push_back(obsy0[i] - a);
		obsv.push_back(obsvx0[i]);
		obstheta.push_back(0.0);
	}

    t[0] = 0;
    start = clock();
    for (z = 0; z <= N; z++) {
        for (i = 0; i < 1; i++) {
            V1[i] = 0;
            V2[i] = 0;
            viscx[i] = 0;
            viscy[i] = 0;
            viscLx[i] = 0;
            viscLmax[i] = 0;

            // POTENTIALS NUDGING AND VISCOCITY
            for (j = 0; j < obs_n; j++) {
				
				dis[i][j] = distance(x, obsx[j], y, obsy[j], p);
				Vdot[i][j] = Vfunc(dis[i][j], lambda, L, q);
				V1[i] = V1[i] + Vdot[i][j] * (x - obsx[j]) / dis[i][j];
				V2[i] = V2[i] + Vdot[i][j] * (y - obsy[j]) / dis[i][j];

				if (visc != 0.0) { 
				    //FOR viscocity
				    viscx[i] = viscx[i] + kappa(dis[i][j], L, lambda, visc) * (obsv[j] * cos(obstheta[j]) - obsv[i] * cos(obstheta[i]));
				    viscy[i] = viscy[i] + kappa(dis[i][j], L, lambda, visc) * (obsv[j] * sin(obstheta[j]) - obsv[i] * sin(obstheta[i]));
				}
            }
        
            // CONTROLLERS FOR EACH MODEL
            if (model == 1) {
                //NEWTONIAN MODEL
                Li = V1[i] - viscx[i];
                k = m2 + (1 / vdx) * (Li)+((vmax * cos(theta)) / (vdx * (vmax * cos(theta) - vdx))) * func1(-Li, epsilon); // New controller with g_1=g_2=s
                F = -(1 / cos(theta)) * (k * (v * cos(theta) - vdx) + V1[i] - viscx[i]);
                u = -(1 / (vdx + ((double)A / (v * (pow((cos(theta) - cos(phi)), 2)))))) * (m1 * v * sin(theta) + Ufunc(y, a, c, pot) + p * V2[i] + sin(theta) * F - viscy[i]);
			}
            else if (model == 2) {
                // RELATIVISTIC MODEL
                F = 2 * pow((vmax - v), 2) * pow(v, 2) / (vmax * v * cos(theta) - 2 * vdx * v + vdx * vmax) * (viscx[i] - m2 * (v * cos(theta) - vdx) - V1[i]);
                u = v * (viscy[i] - m1 * v * sin(theta) - Ufunc(y, a, c, pot) - vmax * sin(theta) * F / (2 * pow((vmax - v), 2) * v) - p * V2[i]) / (A / pow(cos(theta) - cos(phi), 2) + vdx / (vmax - v));
			}
            delta = atan(u * sigma / v);

            if (Fmax <= F) {
                Fmax = F;
            }

            //SAMPLED DATA MODEL
            if (delta == 0.0) {
                x = x + (v * cos(theta)) * T + F * cos(theta) * pow(T, 2) / 2;
                y = y + (v * sin(theta)) * T + F * sin(theta) * pow(T, 2) / 2;
                // theta = theta;
                v = v + T * F;
            }
            else {
                x = x + sigma / tan(delta) * (sin(theta + v * tan(delta) * T / sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma)) - sin(theta));
                y = y + sigma / tan(delta) * (cos(theta) - cos(theta + v * tan(delta) * T / sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma)));
                theta = theta + v * tan(delta) * T / (double)sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma);
                v = v + T * F;
            }

            if ((fabs(y) >= a) || (fabs(theta) > phi)) {
                flag_fail = 1;
                printf(RED "FAIL:\t" RESET "Collision occured\t % f\t % f\n", theta, y);
                break;
            }

			// calculate obstacle trajectories
			for (int o = 0; o < obs_n; o++) {
                F = 0;
                obsx[o] = obsx[o] + (obsv[o] * cos(obstheta[o])) * T + F * cos(obstheta[o]) * pow(T, 2) / 2;
                obsy[o] = obsy[o] + (obsv[o] * sin(obstheta[o])) * T + F * sin(obstheta[o]) * pow(T, 2) / 2;
                obstheta[o] = obstheta[o];
                obsv[o] = obsv[o] + T * F;
            }
        }
        if (flag_fail == 1) break;
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    if (flag_fail == 0) {
        printf("time: %.10f \n", cpu_time_used);
		printf(GRN "SUCCESS:" RESET "%d\t%lf\t%lf\t%lf\t%lf\n", z, v, theta, y, x);
    }
}

