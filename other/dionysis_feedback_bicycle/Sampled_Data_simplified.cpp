#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include <vector>

using namespace std;

#define MAX(a,b) ((a > b)?(a):(b))

/* Colors */
#define RED   "\x1B[1;31m"
#define GRN   "\x1B[1;32m"
#define RESET "\x1B[0m"


/* Parameters */
#define a 5.0           // perpendicular of the road in meters (half road width)
#define sigma 5         // Vehicle length in meters
#define phi 0.25        // Angle max
// #define p 5.11251       // 1/(3*pow(tan(phi),2))
#define p 15.0       // 1/(3*pow(tan(phi),2))
#define L 5.5       // sigma*fmax(2*sqrt(p)*sin(phi),sqrt(1+(p-1)*pow(sin(phi),2)))


#define epsilon 0.2
#define lambda 25       // radius of communication
#define c  1.5          // 3.86866 //1.5
#define vmax 35         // Maximum velocity of vehicles

#define vstar 30        // set point velocity
#define A 1
int n;             // number of vehicles(+1)


#define pot 0

#define model 1         // 1. Newtonian -- 2. Relativistic

#if (model == 1)
    // MODEL NEWTONIAN

    #define m1 0.5
    // #define m2   1/vmax
    #define m2   0.5
    double qu = 1;
    double q = 0.003;
#elif (model == 2)
    // MODEL RELATIVISTIC

    #define m1 0.3
    #define m2 1/pow(vmax,2)
    double qu =  1;
    double q = 0.000003;
#endif

// VISOSCITY      
//double visc = 0.05/pow(vmax,2);      
// double visc =  1 / vmax;
double visc = 0.0;   // 0 if we dont want to take other vehicle speeds into account      

int tmax = 20; //  Time ending point         
double T = 0.0015; // Step-size  
int flag_fail = 0;

int every = 1000;


int counT = 1;
double Fmax = 0;

double distance(double xi, double xj, double yi, double yj, double var);
double Vfunc(double var1, int var2, double var3, double var4);
double sum(double(*), int var2);
double func1(double var2, double var1);
double Ufunc(double var1, double var2, double var3, int var4);
double kappa(double var1, double var2, double var3, double var4);
int mod(int v1, int v2);

struct Vehicle {
    vector<double> vdx, x, y, theta, vx;
};

struct Parameters{
    int numsteps;
    double step;
    vector<double> vdx, x, y, theta, vx;
}P;

void read() {

	char buf[1024];
	// const char* line;
	double dval;
	int ival, numVeh = 0;


	FILE* input;
	input = fopen("input.txt", "r");

	while (fgets(buf, sizeof(buf), input)) {
		// fprintf(stderr, "parsing: %s\n", buf);
		if (buf[0] == '\n')
			break;

		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			P.numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			P.step = dval;
		
		/* read initial states */
		if (sscanf(buf, "\"vdx(%d)\":%lf", &ival, &dval) == 2)
			P.vdx.push_back(dval);
        if (sscanf(buf, "\"x(%d)\":%lf", &ival, &dval) == 2)
			P.x.push_back(dval);
        if (sscanf(buf, "\"y(%d)\":%lf", &ival, &dval) == 2)
			P.y.push_back(dval);
        if (sscanf(buf, "\"vx(%d)\":%lf", &ival, &dval) == 2)
			P.vx.push_back(dval);
        if (sscanf(buf, "\"theta(%d)\":%lf", &ival, &dval) == 2)
			P.theta.push_back(dval);
		
		memset(buf, 0, sizeof(buf));
	}

	n = P.x.size();
}

void exportSimData(struct Vehicle *veh, int numsteps, double step, int numVeh) {
	FILE* fd2 = NULL;
	if ((fd2 = fopen("../../sim/data/sim.js", "w")) == NULL) { printf("\nCouldn't export data file \n"); }

	double max_pos = 0.0;
	fprintf(fd2, "sim = {\n");
	fprintf(fd2, "x: [\n");
	
	// longitudinal position
	double pos_x;
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			pos_x = veh[idVeh].x[k];
			fprintf(fd2, "%.4f,", pos_x);
			if (pos_x > max_pos)
				max_pos = pos_x;
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	// lateral position
	double pos_y;
	fputs("y: [\n", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			pos_y = veh[idVeh].y[k] + a;
			fprintf(fd2, "%.4f,", pos_y);
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	// longitudinal speed
	double spd_x;
	fputs("vx: [\n", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			spd_x = veh[idVeh].vx[k];
			fprintf(fd2, "%.4f,", spd_x);
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("vdx: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%.4f,", P.vdx[idVeh]);
	fputs("],\n", fd2);

	fputs("id: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%d,", idVeh + 1);
	fputs("],\n", fd2);

	fputs("Cx: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < numVeh; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", numVeh);
	fprintf(fd2, "k:%d,\n", numsteps);
	fprintf(fd2, "roadlength:%f,\n", max_pos);
	fprintf(fd2, "roadwidth:%f,\n", 10.0);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", step);
	// fprintf(fd2, "safety1:%f,\n", this->safety_x);
	// fprintf(fd2, "safety2:%f,\n", this->safety_y);
	fprintf(fd2, "roadbound:%f,\n", 1.2);

	fprintf(fd2, "};\n\n");

	fclose(fd2);
}

int main(void) {

    Vehicle *veh;
    read();
    veh = (struct Vehicle*)calloc(n, sizeof(struct Vehicle));
    double **init;
    init = (double**)calloc(4, sizeof(double*));
    for (int i = 0; i < 4; i++)
        init[i] = (double*)calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        // fscanf(myFile, "%lf", &init[s][i]);
        init[0][i] = P.x[i];
        init[1][i] = P.y[i] - a;
        init[2][i] = P.theta[i];
        init[3][i] = P.vx[i];

        // asdasd
        veh[i].x.push_back(P.x[i]);
        veh[i].y.push_back(P.y[i] - a);
        veh[i].theta.push_back(P.theta[i]);
        veh[i].vx.push_back(P.vx[i]);
    }

    double *mindT, *t, **x, **y, **v, **theta, **xf, **yf, **thetaf, **vf;

    clock_t start, end;
    double cpu_time_used;

    double delta;
    int  i, j, z;
    int N = (int)(tmax / T);
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

    x = (double**)calloc(1, sizeof(double*));
    y = (double**)calloc(1, sizeof(double*));
    v = (double**)calloc(1, sizeof(double*));
    theta = (double**)calloc(1, sizeof(double*));
    mindT = (double*)calloc(1, sizeof(double));
    xf = (double**)calloc(N2 + 1, sizeof(double*));
    yf = (double**)calloc(N2 + 1, sizeof(double*));
    thetaf = (double**)calloc(N2 + 1, sizeof(double*));
    vf = (double**)calloc(N2 + 1, sizeof(double*));
    for (i = 0; i < 1; i++) {
        x[i] = (double*)calloc(n, sizeof(double));
        y[i] = (double*)calloc(n, sizeof(double));
        v[i] = (double*)calloc(n, sizeof(double));
        theta[i] = (double*)calloc(n, sizeof(double));
    }
    for (i = 0; i < N2 + 1; i++) {
        xf[i] = (double*)calloc(n, sizeof(double));
        yf[i] = (double*)calloc(n, sizeof(double));
        vf[i] = (double*)calloc(n, sizeof(double));
        thetaf[i] = (double*)calloc(n, sizeof(double));
    }
    t = (double*)calloc(N2 + 1, sizeof(double));

    for (j = 0; j < n; j++) {
        x[0][j] = init[0][j];
        xf[0][j] = init[0][j];
        y[0][j] = init[1][j];
        yf[0][j] = init[1][j];
        theta[0][j] = init[2][j];
        thetaf[0][j] = init[2][j];
        v[0][j] = init[3][j];
        vf[0][j] = init[3][j];
    }

    t[0] = 0;
    start = clock();
    for (z = 0; z <= N; z++) {
        for (i = 0; i <= n - 1; i++) {
            
            if (i != 0) {
                F = 0;
                x[0][i] = x[0][i] + (v[0][i] * cos(theta[0][i])) * T + F * cos(theta[0][i]) * pow(T, 2) / 2;
                y[0][i] = y[0][i] + (v[0][i] * sin(theta[0][i])) * T + F * sin(theta[0][i]) * pow(T, 2) / 2;
                theta[0][i] = theta[0][i];
                v[0][i] = v[0][i] + T * F;

                veh[i].x.push_back(x[0][i]); veh[i].y.push_back(y[0][i]);
                veh[i].theta.push_back(theta[0][i]); veh[i].vx.push_back(v[0][i]);

                continue;
            }

            V1[i] = 0;
            V2[i] = 0;
            viscx[i] = 0;
            viscy[i] = 0;
            viscLx[i] = 0;
            viscLmax[i] = 0;

            // POTENTIALS NUDGING AND VISCOCITY
            for (j = 0; j <= n - 1; j++) {
                if (i != j) {
                    dis[i][j] = distance(x[0][i], x[0][j], y[0][i], y[0][j], p);
                    Vdot[i][j] = Vfunc(dis[i][j], lambda, L, q);
                    V1[i] = V1[i] + Vdot[i][j] * (x[0][i] - x[0][j]) / dis[i][j];
                    V2[i] = V2[i] + Vdot[i][j] * (y[0][i] - y[0][j]) / dis[i][j];

                    if (visc != 0.0) { 
                        //FOR viscocity
                        viscx[i] = viscx[i] + kappa(dis[i][j], L, lambda, visc) * (v[0][j] * cos(theta[0][j]) - v[0][i] * cos(theta[0][i]));
                        viscy[i] = viscy[i] + kappa(dis[i][j], L, lambda, visc) * (v[0][j] * sin(theta[0][j]) - v[0][i] * sin(theta[0][i]));
                    }
                }
            }
        
            // CONTROLLERS FOR EACH MODEL
            #if (model == 1)
                //NEWTONIAN MODEL
                Li = V1[i] - viscx[i];
                k = m2 + (1 / P.vdx[i]) * (Li)+(((double)vmax * cos(theta[0][i])) / (P.vdx[i] * ((double)vmax * cos(theta[0][i]) - P.vdx[i]))) * func1(-Li, epsilon); // New controller with g_1=g_2=s
                F = -(1 / cos(theta[0][i])) * (k * (v[0][i] * cos(theta[0][i]) - P.vdx[i]) + V1[i] - viscx[i]);
                u = -(1 / (P.vdx[i] + ((double)A / (v[0][i] * (pow((cos(theta[0][i]) - cos(phi)), 2)))))) * (m1 * v[0][i] * sin(theta[0][i]) + Ufunc(y[0][i], a, c, pot) + p * V2[i] + sin(theta[0][i]) * F - viscy[i]);
            #else
                // RELATIVISTIC MODEL
                F = 2 * pow((vmax - v[0][i]), 2) * pow(v[0][i], 2) / (vmax * v[0][i] * cos(theta[0][i]) - 2 * P.vdx[i] * v[0][i] + P.vdx[i] * vmax) * (viscx[i] - m2 * (v[0][i] * cos(theta[0][i]) - P.vdx[i]) - V1[i]);
                u = v[0][i] * (viscy[i] - m1 * v[0][i] * sin(theta[0][i]) - Ufunc(y[0][i], a, c, pot) - vmax * sin(theta[0][i]) * F / (2 * pow((vmax - v[0][i]), 2) * v[0][i]) - p * V2[i]) / (A / pow(cos(theta[0][i]) - cos(phi), 2) + P.vdx[i] / (vmax - v[0][i]));
            #endif
            delta = atan(u * sigma / v[0][i]);

            if (Fmax <= F) {
                Fmax = F;
            }

            //SAMPLED DATA MODEL
            if (delta == 0.0) {
                x[0][i] = x[0][i] + (v[0][i] * cos(theta[0][i])) * T + F * cos(theta[0][i]) * pow(T, 2) / 2;
                y[0][i] = y[0][i] + (v[0][i] * sin(theta[0][i])) * T + F * sin(theta[0][i]) * pow(T, 2) / 2;
                theta[0][i] = theta[0][i];
                v[0][i] = v[0][i] + T * F;
            }
            else {
                x[0][i] = x[0][i] + sigma / tan(delta) * (sin(theta[0][i] + v[0][i] * tan(delta) * T / sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma)) - sin(theta[0][i]));
                y[0][i] = y[0][i] + sigma / tan(delta) * (cos(theta[0][i]) - cos(theta[0][i] + v[0][i] * tan(delta) * T / sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma)));
                theta[0][i] = theta[0][i] + v[0][i] * tan(delta) * T / (double)sigma + F * tan(delta) * pow(T, 2) / (2 * (double)sigma);
                v[0][i] = v[0][i] + T * F;
            }

            if ((fabs(y[0][i]) >= a) || (fabs(theta[0][i]) > phi)) {
                flag_fail = 1;
                printf(RED "FAIL:\t" RESET);
                printf("Some collision occured\t % f\t % f\n", theta[0][i], y[0][i]);
                break;
            }

            // keep data for sim
            veh[i].x.push_back(x[0][i]); veh[i].y.push_back(y[0][i]);
            veh[i].theta.push_back(theta[0][i]); veh[i].vx.push_back(v[0][i]);
        }
        if (flag_fail == 1) {
            printf("%d", z);

            break;
        }
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    exportSimData(veh, N, T, n);
    
    if (flag_fail == 0) {
        printf("time: %.10f \n", cpu_time_used);
        if (vf[N2][1] != 0) {
            printf("%d\t%lf\t%lf\t%lf\t%lf\n", z, vf[N2][1], thetaf[N2][1], yf[N2][1], xf[N2][1]);
        }
        else {
            printf(GRN "SUCCESS:" RESET "%d\t%lf\t%lf\t%lf\t%lf\n", z, v[0][1], theta[0][1], y[0][1], x[0][1]);
        }
    }

    return 0;
}


double sum(double(*var1), int var2) {
    int var4;
    double result = 0;
    for (var4 = 0; var4 <= var2; var4++)
    {
        result += var1[var4];
    }
    return result;
}

double distance(double xi, double xj, double yi, double yj, double var) {
    double result;
    //result = sqrt((pow((xi - xj),2) + var*pow((yi - yj),2)));
    result = sqrt((xi - xj) * (xi - xj) + p * (yi - yj) * (yi - yj));
    return result;
}

double Vfunc(double var1, int var2, double var3, double var4) {
    double result;
    if (var1 > var2) {
        result = 0;
    }
    else {
        result = (double)(var4 * (-(double)3 * pow(((double)var2 - var1), 2) * (var1 - var3) - pow(((double)var2 - var1), 3))) / (pow((var1 - var3), 2));
    }
    
    return result;
}

double func1(double var2, double var1) {
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

double Ufunc(double var1, double var2, double var3, int var4) {
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

double kappa(double var1, double var2, double var3, double var4) {
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

int mod(int v1, int v2) {
    int ret = v1 % v2;
    if (ret < 0)
        ret += v2;
    return ret;
}