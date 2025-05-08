#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>
#include <string.h>

#include <vector>

using namespace std;

#define CONTROLLER_H
#include "main.h"

#define BUS_RATIO 0
double LANEWIDTH;
#define VD_VARIANCE 2.0
#define XBUFFER 2.0

#define SAMPLE_UNIFORM(min, max) ((double)min + ((double)random()/RAND_MAX)*(max - min))
#define MAX(a, b) (((a) > (b))?(a):(b))
#define MIN(a, b) (((a) <= (b))?(a):(b))

int deiktis3;
static int pl = 0;
static double pl2 = 0;
static double pl3 = 0;

#define EPSILON_SIZE 6

#define P_CR 360 

#define BUFFER_SIZE_CHAR_P_ROW 30*(3*EPSILON_SIZE) 
#define D1 5000
#define D2 200 

struct Vehicle {
	double x, y, vx, vy, vdx, vdy, w, l;
}veh = {0};

struct Obstacles {
	double x, y, vx, vy, vdx, vdy, w, l;
}obs = {0};

struct {
	int numsteps;
	double step;
	double vdx, vdy, safety_x, safety_y, C;
	double UXMIN, UXMAX, UYMIN, UYMAX;
	double VXMIN, VXMAX;
	double x0, y0, vx0, vy0;

	vector<double> obs_x0, obs_y0, obs_vx0, obs_vy0;
	int obs_n;

} P = { 0 };

// longitudinal target speed
static double ftsx(sim_t* sim, double vx, double vd) {  
	// trying to get the desired speed

	return 1 * (erfc(0.2 * (vx - vd)) - 1);
}

//lateral target speed
static double ftsy(sim_t* sim, double vy, double vdy) { 
	// trying to get the desired speed
	
	//return sim->ftsy_hi*(2 / (1 + exp(-sim->ftsy_zeta*(sim->vdy_effective[i] - sim->vy[i][t]))) - 1);
	return 1 * (erfc(0.5 * (vy - vdy)) - 1);
}

//G and H (below are used for estimating the lateral aura.
static double G(double t, double warmup, double length, double cooldown) {
    return fmax(0, fmin(fmin(1. , t/warmup), 1. - (t-warmup-length)/cooldown));
}

static double H(double t, double warmup, double length, double cooldown) {
    return G(t+warmup+length, warmup, 2*length, cooldown);
}
 
//calculate the forces between our vehicle and one of its neighbours
static long double fca_mag(sim_t* sim, NumericalID edge_id, double dx, double dy, double vxi, double vxj, double approaching_vx, double approaching_vy, double len, double wid, NumericalID veh_i, NumericalID veh_j, double* ca) {
	
	double magx = 0;
	double Rdesired;
	double Rmin;
	double dx3;
	double dv2;
	int pointer = 0;
	int pointer2 = 0;

	double Dmin = sim->D_min; //  0.5; // 1 input file
	double Dmax = sim->D_max;// 2.5; // 3 input file diam  2.5
	double Dmax2 = sim->D_max2;// 2.7;  //2.7
	double yi = get_position_y(veh_i);
	double xi = get_position_x(veh_i);
	double Dy = wid;
	double codiag = 100;
	double Coef = 0.4;
	double ds = 1.2;

	double Ly = approaching_vy * sim->time_gap_y + 0.65;// Ly from page 12 3d line 

	double magy = H(dy, Ly, Dy + 0.65, Ly);

	if (magy == 0) {
		pointer = 0;    // if the neighbor is not an obstacle
	}
	else {
		pointer = 1;    // if the neighbor is an obstacle (is in the lateral aura)
	}

	double dx2 = dx;

	// if j approaches i
	if (dx2 < 0) {
		pointer2 = 1;
		dx2 = -dx2;

		Rdesired = ds + len + vxj * sim->time_gap_x;
		Rmin = ds + len + Coef * vxj * sim->time_gap_x;// 0.5*Rdesired;

		dv2 = vxj - vxi;
	}
	else {
		// if i approaches j  normal situation
		pointer2 = 0;

		Rdesired = ds + len + vxj * sim->time_gap_x;
		Rmin = ds + len + Coef * vxj * sim->time_gap_x;// 0.5*Rdesired;

		dv2 = vxj - vxi;
	}

	
	double Rsw = 0, Rmn = 0, sp = 2;
	//////////////////////////////////////////////////////////////////////////////////////
	double curve3 = Rdesired + (pow((dv2), 2) / (2 * Dmin));
	double curve1 = Rdesired + (pow((dv2), 2) / (2 * Dmax));
	double curve2 = Rmin + (pow((dv2), 2) / (2 * Dmax2));
	double curve4 = Rmin - (Rmin / sp) * (dv2);

	deiktis = 0;
	deiktis3 = 0;
	/////////////////////////////////////////////////////////////////////////////////////	
	*ca = 1;
	/////////////////////////////////////////////////

	if (pointer == 1) {
		if (dv2 >= 0) {
			magx = 0;
			deiktis = 1;
		}

		if ((dv2 < 0) && (dx2 > curve2)) {
			magx = pow((dv2), 2) / (2 * (dx2 - Rmin));

			if (pointer2 == 0) {
				deiktis = 3; // between the min nad max decelerations that leads to Rdesired
			}
			else {
				deiktis3 = 3;
			}
		}

		if ((dv2 < 0) && (dx2 <= curve2)) {
			magx = 2.5;

			if (pointer2 == 0) {
				deiktis = 5;  // emergency deceleration
			}
			else {
				deiktis3 = 5;
			}
		}
		if ((dv2 >= 0) && (dv2 < sp) && (dx2 <= curve4)) {

			magx = 2.5;

			if (pointer2 == 0) {
				deiktis = 55;  // emergency deceleration
			}
			else {
				deiktis3 = 55;
			}
		}
	}
	else {
		magx = 0;
		deiktis = 11;
		deiktis3 = 11;
	}
	
	if (pointer2 == 1)  {
		//if j approaches i 
		magx = -magx;
	}

	return magx * magy;
}

static double fca(sim_t *sim, NumericalID edge_id, NumericalID veh_i, NumericalID veh_j, double dx_i_to_j, double dy_i_to_j, double *fcax, double *fcay, double* caa) {
    double mag, ang, ca;
	double li = get_veh_length(veh_i), lj = get_veh_length(veh_j);
	double wi = get_veh_width(veh_i), wj = get_veh_width(veh_j);
    double len = 0.5*(li + lj);
    double wid = 0.5*(wi + wj);

	*fcax = *fcay = *caa = 0;

	double x, x0, y, y0, approaching_speed_x, approaching_speed_y;
    x0 = get_position_x(veh_j);  
    y0 = get_position_y(veh_j);
    x = get_position_x(veh_i); 
    y = get_position_y(veh_i); 
	 
	double vxi = get_speed_x(veh_i), vxj = get_speed_x(veh_j);
	double vyi = get_speed_y(veh_i), vyj = get_speed_y(veh_j);

		
	approaching_speed_x = MAX(0, vxi - (1. - sim->safety_level_x) * vxj);
	approaching_speed_y = MAX(0, ((vyi - vyj) * ((y0 - y) / fabs(y0 - y))));

    mag = fca_mag(sim, edge_id,
                dx_i_to_j, dy_i_to_j, 
                vxi,vxj,
                approaching_speed_x, approaching_speed_y, len, wid, veh_i, veh_j, & ca);
    
    ang = atan2(dy_i_to_j, fabs(dx_i_to_j));
    *fcax = -cos(ang) * mag;
    *fcay = -sin(ang) * mag;
	*caa = ca;

    return mag;
}

int cmpnbors(const void *n1ptr, const void *n2ptr) {
    nbor_t *n1 = (nbor_t*)n1ptr;
    nbor_t *n2 = (nbor_t*)n2ptr;

    if (n1->mag == n2->mag)
        return 0;
    if (n1->mag < n2->mag)
        return -1;
    
    return +1;
}

void determine_forces(sim_t* sim, struct Vehicle v, struct Obstacles *obs) { 
   
	static double taccsm[30000];
	
	int j;
    double fcax = 0, fcay = 0;
    double vdy2 = 0;
	
	int obs_down = 0, obs_up = 0;
	for (int i = 0; i < P.obs_n; i++) {
		if (obs[i].x >= v.x)
			obs_down++;
		else
			obs_up++;
	}


    // int max_n = n;
    // if (n_neighbors_front != -1)
	// 	max_n = n_neighbors_front;

	// int max_n_back = n;
	// if (n_neighbors_back != -1)
	// 	max_n_back = n_neighbors_back;
	
    nbor_t* nbors = NULL; 
	if (obs_up > 0)
		nbors = (nbor_t*)malloc(obs_up * sizeof(nbor_t));

	nbor_t* nbors_nudge = NULL;
	if (obs_down > 0)
		nbors_nudge = (nbor_t*)malloc(obs_down * sizeof(nbor_t));
	
	
	// infrastructure parameters
	double roadwid_meters = 10.0;

	// CHANGED
	/* ~~~~~~~~~ set desired speed and vdy2(what is this???) -- I commented all and set vdy2 = 0 ~~~~~~~~~~~~~~~~~~~ */
	vdy2 = 0;	// wtf is this
	/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
	
	double fxi = 0, fyi = 0, caa = 1;
	
    /* obstacle-related forces */
    double fcax_sum = 0, fcay_sum =0;
    int nbors_added = 0, nbors_added_nudge = 0;
	
    
    int deiktis = 0, deiktis2 = 0;
    int correction = 0;
	double Ma = 0;
	
	// one for from i until front_end, where front_end is n-1, and when it reaches n-1 it resets to zero. we break when either the ifluence_radius_meters 
	// is exceeded or we have relative distance >roadlength/2
	double dx_i_to_j, dy_i_to_j, fcamag;

    NumericalID* downstream_vehs_array;
    int downstream_start_index, downstream_n;
	
	// WTF IS THIS PART
	//Downstream vehicles
    if (n_neighbors_front == -1) { // if n_neighbors==-1, then we do not use the functionality of front neighbors, so for downstream vehicles we retain the old functionality
        downstream_start_index = i + 1;
        downstream_vehs_array = vehs_array;
        downstream_n = n;
    }
    else { // use the front neighbors array, and start from 0
        downstream_start_index = 0;
        downstream_vehs_array = front_neighbors;
        downstream_n = n_neighbors_front;
    }

	//Upstream vehicles
	NumericalID* upstream_vehs_array;
	int upstream_start_index, upstream_n;
	if (n_neighbors_back == -1) {//if n_neighbors_back==-1, then we do not use the functionality of back neighbors, so for upstream vehicles we retain the old functionality
		upstream_start_index = i - 1;
		upstream_vehs_array = vehs_array;
		upstream_n = n;
	}
	else { // use the front neighbors array, and start from 0
		upstream_start_index = 0;
		upstream_vehs_array = back_neighbors;
		upstream_n = n_neighbors_back;
	}

    for (j = downstream_start_index; ; j ++) {
        //printf("%d, ", j);
        if (n_neighbors_front == -1) {
            if (j > downstream_n - 1)
                j = 0;

            if (j == i)
                break;
        }
        else {
            if (j >= n_neighbors_front)
                break;
        }
        
        dx_i_to_j = get_relative_distance_x(veh_id, downstream_vehs_array[j]);
        dy_i_to_j = get_relative_distance_y(veh_id, downstream_vehs_array[j]);
	
        fcamag = 0;
        
        if (n_neighbors_front == -1) {
            if ((dx_i_to_j) > MIN(sim->influence_radius_meters, 1000) || dx_i_to_j < 0)
                break;
        }		

        fcamag = fca(sim, edge_id, veh_id, downstream_vehs_array[j], dx_i_to_j, dy_i_to_j, &fcax, &fcay, &caa);

        if ((fcax <= 0)) {
            nbors[nbors_added].j = j;
            nbors[nbors_added].mag =  -fcamag;
            nbors[nbors_added].fcax = fcax;
            nbors[nbors_added].fcay = fcay;
			nbors[nbors_added].deiktis4 = deiktis;

            nbors_added++;
        }   
	}
      
		
	for (j = upstream_start_index; ; j = (j - 1) * (n_neighbors_back == -1) + (j + 1) * (n_neighbors_back != -1)) {
		if (n_neighbors_back == -1) {
			if (j < 0)
				j = upstream_n - 1;
	
			if (j == i)
				break;
		}
		else {
			if (j >= n_neighbors_back)
				break;
		}

		dx_i_to_j = get_relative_distance_x(veh_id, upstream_vehs_array[j]);
		dy_i_to_j = get_relative_distance_y(veh_id, upstream_vehs_array[j]);
		fcamag = 0;
		
		if ((-dx_i_to_j) > MIN(sim->influence_radius_meters, 1000) || dx_i_to_j > 0)
			break;

		fcamag = fca(sim, edge_id, upstream_vehs_array[j], veh_id, dx_i_to_j, -dy_i_to_j, &fcax, &fcay , & caa);
		
		if (fcax > 0) {
			nbors_nudge[nbors_added_nudge].j = j;
			nbors_nudge[nbors_added_nudge].mag =  fcamag;
			nbors_nudge[nbors_added_nudge].fcax = fcax;
			nbors_nudge[nbors_added_nudge].fcay = fcay;
			nbors_nudge[nbors_added_nudge].deiktis5 = deiktis3;
			nbors_added_nudge++;
		}	
	}

    // sort neighbors and compute fcax_sum from top sim->nbors
    if (sim->fca_nbors && nbors_added > 0) {
        int x;
        qsort(nbors, nbors_added, sizeof(nbor_t), cmpnbors);

		deiktis2 = nbors[0].deiktis4;
		
		Ma = -nbors[0].mag;
		if (deiktis2 <= 0)
			deiktis2 = 1;
		
        fcax_sum = fcay_sum = 0;
		double y_nbors; double x_nbors;
		for (x=0; x < MIN(sim->fca_nbors, nbors_added); x++) {
			int jj = nbors[x].j;
			y_nbors = get_position_y(downstream_vehs_array[nbors[x].j]);
			fcax_sum += nbors[x].fcax;

			if ((strcmp(get_edge_name(edge_id), "segment5") == 0) || (strcmp(get_edge_name(edge_id), "segment3_opp") == 0)) {
				if ((y_nbors <= 3 - wi / 2) && (xi > 100))
					fcay_sum += nbors[x].fcay;
			
				if (y_nbors > 3 - wi / 2)
					fcay_sum +=  nbors[x].fcay;
			}
			else {
				fcay_sum += nbors[x].fcay;
			}

		}
	}

    if (sim->fca_nbors_nudge && nbors_added_nudge > 0) {
        int x;
		qsort(nbors_nudge, nbors_added_nudge, sizeof(nbor_t), cmpnbors);
		
        for (x=0; x < MIN(sim->fca_nbors_nudge, nbors_added_nudge); x++) {
				
			double y_nbors_nudging = get_position_y(vehs_array[j]);

			if (Ma <= 2) {
				fcax_sum += nbors_nudge[x].fcax * sim->fwd_force_max_x;

				if ((strcmp(get_edge_name(edge_id), "segment5") == 0) || (strcmp(get_edge_name(edge_id), "segment3_opp") == 0)) {

					if ((y_nbors_nudging <= 3 - wi / 2) && (xi > 100)) {
						fcay_sum += nbors_nudge[x].fcay * sim->fwd_force_max_y; 
					}
					if (y_nbors_nudging > 3 - wi / 2)
					{
						fcay_sum += 1 * nbors_nudge[x].fcay * sim->fwd_force_max_y;
					}
				}
				else {
					fcay_sum += nbors_nudge[x].fcay * sim->fwd_force_max_y;
				}
			}

		}
	}
	
	double vdii2;

	fxi += sim->coeff_fcax * fcax_sum;
	fyi += sim->coeff_fcay * fcay_sum;

	double vxi = get_speed_x(veh_id), vyi = get_speed_y(veh_id), vdii = get_desired_speed(veh_id);
    /* target-speed related forces */

	double vdi = vdii;//  MIN(vxi + 7, vdii);

    if ((deiktis2 == 1) || (deiktis2 == 11) || (deiktis2 == 0)) //speed control teritory
        fxi += ftsx(sim, vxi, vdii);
	

    fyi += ftsy(sim, vyi, vdy2);


    if ((fxi > taccsm[veh_id]) && (taccsm[veh_id]>=0))
        taccsm[veh_id] = 0.4 * fxi + 0.6 * taccsm[veh_id];
    else
        taccsm[veh_id] = 0.6 * fxi + 0.4 * taccsm[veh_id];
    

    // if (((xi <= 150) && (strcmp(get_edge_name(edge_id), "warm_up") == 0)) || ((xi <= 150) && (strcmp(get_edge_name(edge_id), "warm_up_opp") == 0))) {
    //     taccsm[veh_id] = 1*ftsx(sim, vxi, vdi); // diam +
    // }



    *fx = taccsm[veh_id];
    *fy = fyi;
    *vdy = vdy2;

    free(nbors);
    free(nbors_nudge);
}

void regulate_forces(sim_t* sim, NumericalID edge_id, NumericalID veh_id, double vdy, double* fx, double* fy, FILE* M4) {
	double fxi = *fx, fyi = *fy;

	double vx = get_speed_x(veh_id), vd = get_desired_speed(veh_id), vy = get_speed_y(veh_id);
	double yi = get_position_y(veh_id), wi = get_veh_width(veh_id), li = get_veh_length(veh_id);
	double xi = get_position_x(veh_id);
	double roadwid_meters = get_edge_width(edge_id);
	double T = get_time_step_length();


	// the acceleration lane's vehicle breaks to avoid an accident
	 if ((strcmp(get_edge_name(edge_id), "segment5") == 0) || (strcmp(get_edge_name(edge_id), "segment3_opp") == 0)) {
		if (yi < 3 + wi/2  && xi>50) {
			if (vdy > 0) {
				if ((vy + fyi * T)<= 0) {
					double yy = (vy + fyi * T)/vdy;
					if (yy < -1) {
						yy = -1;
					}
					if (fxi > 0) {
						fxi = fxi -1.5+ 1.5* (yy - (-1)); ///min value -1.5 at -1 max value 0.0  at 0.00					
					}
				}
			}
		}
	}

	/* [umin, umax] ranges */
	fxi = MIN(fxi, sim->uxmax_hard);
	fxi = MAX(fxi, sim->uxmin_hard);
	
	/* non-negative speed */
	fxi = MAX(fxi, -vx / T);

	/* non-excessive speed */
	fxi = MIN(fxi, (1.2*vd - vx) / T);
	
	/* keeping vehicles within the road */
	double ydn_next = 0;
	double yup_next = roadwid_meters;
	double vy_bd = 0;
	double vy_bu = 0;
	double ydn = 0; //0
	double yup = roadwid_meters;
	double shibu = 0;
	double acd = 0;
	double acu = 0;
	double shibd = 0;
	double k1 = 0.4;
	double k2 = 1.22;
	double vx_next = vx + T * fxi;
	double x_next = xi + vx*T+0.5* pow(T,2) * fxi;
	
	double distance_from_left_boundary = 0, distance_from_right_boundary = 0, left_boundary_speed = 0, right_boundary_speed = 0;
	get_distance_to_road_boundaries_at(veh_id, 0, 0, &distance_from_left_boundary, &distance_from_right_boundary, &left_boundary_speed, &right_boundary_speed, NULL);
	
	if( (strcmp(get_vehicle_name(veh_id), "highway4.783") == 0) || (strcmp(get_vehicle_name(veh_id), "highway6.231") == 0) || (strcmp(get_vehicle_name(veh_id), "highway7.186") == 0)
		|| (strcmp(get_vehicle_name(veh_id), "highway8.111") == 0) || (strcmp(get_vehicle_name(veh_id), "highway9.27") == 0) || (strcmp(get_vehicle_name(veh_id), "highway10.145") == 0)
		|| (strcmp(get_vehicle_name(veh_id), "highway11.39") == 0)) {
		double left_boundary_global_position = 0, right_boundary_global_position = 0;
		get_global_position_of_road_boundaries_at(veh_id, 0, &left_boundary_global_position, &right_boundary_global_position);
		fprintf(M4, "%lld\t\t%f\t\n", veh_id, left_boundary_global_position);        fflush(M4);
		
	}

	double V_max;
	fyi = MIN(fyi, +(0.05 * vx - vy) / T); 
	fyi = MAX(fyi, -(0.05 * vx + vy) / T);
	V_max = 1;
	 
	// maximum lateral speed bound//    
	if (vy >= 0) {  //last change
		fyi = MIN(fyi, (V_max - vy) /T); //2  the minimum allowed lateral speed  �������� �� �� ���� �� �������
		fyi = MAX(fyi, (-V_max - vy) /T);
	}
	else { 
		fyi = MIN(fyi, (V_max - vy) /T); 
		fyi = MAX(fyi, (-V_max - vy) /T);

	}
  

	if (fyi < -k1 * (distance_from_right_boundary - wi / 2) - k2 * (vy - right_boundary_speed)) {
		fyi = MAX(fyi, -k1 * (distance_from_right_boundary - wi / 2) - k2 * (vy - right_boundary_speed));
	}

	if (fyi > -k1 * (-distance_from_left_boundary + wi / 2) - k2 * (vy - left_boundary_speed)) {
		fyi = MIN(fyi, -k1 * (-distance_from_left_boundary + wi / 2) - k2 * (vy - left_boundary_speed));
	}

    fyi = MIN(fyi, +sim->uymax_hard);
    fyi = MAX(fyi, -sim->uymax_hard);
	
	*fx = fxi;
	*fy = fyi;	
}

void determine_controls(sim_t *sim, double *fx, double *fy) {
    *fx = MIN(*fx, sim->uxmax_hard);
    *fy = (*fy >= 0)? MIN(*fy, sim->uymax_hard): MAX(*fy, -sim->uymax_hard);
} 

void read() {
	char buf[1024];
	// const char* line;
	double dval;
	int ival;

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
		if (sscanf(buf, "\"safety_x\":%lf", &dval) == 1)
			P.safety_x = dval;
		if (sscanf(buf, "\"safety_y\":%lf", &dval) == 1)
			P.safety_y = dval;

		/* read initial states */
		if (sscanf(buf, "\"vdx\":%lf", &dval) == 1)
			P.vdx = dval;
		if (sscanf(buf, "\"vdy\":%lf", &dval) == 1)
			P.vdy = dval;
		if (sscanf(buf, "\"x(0)\":%lf", &dval) == 1)
			P.x0 = dval;
		if (sscanf(buf, "\"y(0)\":%lf", &dval) == 1)
			P.y0 = dval;
		if (sscanf(buf, "\"vx(0)\":%lf", &dval) == 1)
			P.vx0 = dval;
		// if (sscanf(buf, "\"ax(0)\":%lf", &dval) == 1)
		// 	P.ax0 = dval;
		// if (sscanf(buf, "\"ay(0)\":%lf", &dval) == 1)
		// 	P.ay0 = dval;
		
		/* obstacles */
		if (sscanf(buf, "\"obst_x(%d,0)\":%lf", &ival, &dval) == 2)
			P.obs_x0.push_back(dval);
		if (sscanf(buf, "\"obst_y(%d,0)\":%lf", &ival, &dval) == 2)
			P.obs_y0.push_back(dval);
		if (sscanf(buf, "\"obst_vx(%d,0)\":%lf", &ival, &dval) == 2)
			P.obs_vx0.push_back(dval);
		if (sscanf(buf, "\"obst_vy(%d,0)\":%lf", &ival, &dval) == 2)
			P.obs_vy0.push_back(dval);

		memset(buf, 0, sizeof(buf));
	}

	if (!(P.obs_x0.size() == P.obs_y0.size() && P.obs_x0.size() == P.obs_vx0.size() && P.obs_x0.size() == P.obs_vy0.size())) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	P.obs_n = P.obs_x0.size();
}

int main() {

	Vehicle veh;
	sim_t *sim;
	Obstacles *obs;

	read();

	return 0;	

	veh.x = 0.0; veh.y = 5.0;
	veh.vx = 20.0; veh.vy = 0.0; veh.vdx = 25.0; veh.vdy = 0.0;
	veh.w = 2.0; veh.l = 5.0;
    
	determine_forces(sim, veh, obs);

    return 0;
}