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
        if (var1 > -var2 + 1 && var1 < var2 - 1) {
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

void Dionysis_Feedback::run(Vehicle veh) {
    double *init;
    double  veh_width = 1.0;        // half veh width, used for the road boundaries (TODO: read it from input)

    init = (double*)calloc(4, sizeof(double));

	init[0] = x0;
	init[1] = y0 - a;
	init[2] = 0.0;	// theta
	init[3] = vx0;
    // store optimal trajectory
    opt_x.push_back(x0); opt_y.push_back(y0); opt_vx.push_back(vx0);

    clock_t start, end;
    double cpu_time_used;

    double delta;
    int  i, j, z;
    N = numsteps; // TODO: change this
 
    double **dis, **Vdot;
    double *V1, *V2, *viscx, *viscy, *viscLx, *viscLmax;

    dis = (double**)calloc(n, sizeof(double*));
    Vdot = (double**)calloc(n, sizeof(double*));
    for (int i = 0; i < n; i++) {
        dis[i] = (double*)calloc(numObs, sizeof(double));
        Vdot[i] = (double*)calloc(numObs, sizeof(double));
    }

    V1 = (double*)calloc(n, sizeof(double));
    V2 = (double*)calloc(n, sizeof(double));
    viscx = (double*)calloc(n, sizeof(double));
    viscy = (double*)calloc(n, sizeof(double));
    viscLx = (double*)calloc(n, sizeof(double));
    viscLmax = (double*)calloc(n, sizeof(double));


    double u = 0, F = 0, k = 0, Li = 0;

	double x = init[0], y = init[1], theta = init[2], v = init[3];
	vector<double> obsx, obsy, obsv, obstheta;

	for (int i = 0; i < numObs; i++) {
		obsx.push_back(obs[i].x[0]);
		obsy.push_back(obs[i].y[0] - a);
		obsv.push_back(obs[i].vx[0]);
		obstheta.push_back(0.0);
	}
    
    start = clock();
    for (z = 0; z <= N; z++) {
        // fprintf(stderr,"%d\n",z);
        for (i = 0; i < n; i++) {
            V1[i] = 0;
            V2[i] = 0;
            viscx[i] = 0;
            viscy[i] = 0;
            viscLx[i] = 0;
            viscLmax[i] = 0;

            // POTENTIALS NUDGING AND VISCOCITY
            for (j = 0; j < numObs; j++) {
				
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
                u = -(1 / (vdx + ((double)A / (v * (pow((cos(theta) - cos(phi)), 2)))))) * (m1 * v * sin(theta) + Ufunc(y, a - veh_width, c, pot) + p * V2[i] + sin(theta) * F - viscy[i]);
			}
            else if (model == 2) {
                // RELATIVISTIC MODEL
                F = 2 * pow((vmax - v), 2) * pow(v, 2) / (vmax * v * cos(theta) - 2 * vdx * v + vdx * vmax) * (viscx[i] - m2 * (v * cos(theta) - vdx) - V1[i]);
                u = v * (viscy[i] - m1 * v * sin(theta) - Ufunc(y, a - veh_width, c, pot) - vmax * sin(theta) * F / (2 * pow((vmax - v), 2) * v) - p * V2[i]) / (A / pow(cos(theta) - cos(phi), 2) + vdx / (vmax - v));
			}
            delta = atan(u * sigma / v);

            if (Fmax <= F) {
                Fmax = F;
            }

            //SAMPLED DATA MODEL
            if (delta == 0.0) {
                x = x + (v * cos(theta)) * step + F * cos(theta) * pow(step, 2) / 2;
                y = y + (v * sin(theta)) * step + F * sin(theta) * pow(step, 2) / 2;
                // theta = theta;
                v = v + step * F;
            }
            else {
                x = x + sigma / tan(delta) * (sin(theta + v * tan(delta) * step / sigma + F * tan(delta) * pow(step, 2) / (2 * (double)sigma)) - sin(theta));
                y = y + sigma / tan(delta) * (cos(theta) - cos(theta + v * tan(delta) * step / sigma + F * tan(delta) * pow(step, 2) / (2 * (double)sigma)));
                theta = theta + v * tan(delta) * step / (double)sigma + F * tan(delta) * pow(step, 2) / (2 * (double)sigma);
                v = v + step * F;
            }

            if ((fabs(y) >= a) || (fabs(theta) > phi)) {
                flag_fail = 1;
                printf(RED "FAIL:\t" RESET "Collision occured\t % f\t % f\n", theta, y);
                break;
            }

            // store optimal trajectory
            opt_x.push_back(x); opt_y.push_back(y + a); opt_vx.push_back(v);

			// calculate obstacle trajectories
			for (int o = 0; o < numObs; o++) {
                F = 0;
                obsx[o] = obsx[o] + (obsv[o] * cos(obstheta[o])) * step + F * cos(obstheta[o]) * pow(step, 2) / 2;
                obsy[o] = obsy[o] + (obsv[o] * sin(obstheta[o])) * step + F * sin(obstheta[o]) * pow(step, 2) / 2;
                obstheta[o] = obstheta[o];
                obsv[o] = obsv[o] + step * F;
            }
        }
        if (flag_fail == 1) break;
    }

    

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    if (flag_fail == 0)
		printf(GRN "SUCCESS: " RESET "z: %d | v: %.4f | theta: %.4f | y: %.4f | x: %.4f | CPU time: %.4f \n", 
            z, v, theta, y, x, cpu_time_used);
}

