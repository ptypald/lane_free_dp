#include <iostream>
#include <vector>
#include <string.h>

// #define DATA_H
#include "data.h"

void Data::read() {

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
			this->numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			this->step = dval;
		if (sscanf(buf, "\"safety_x\":%lf", &dval) == 1)
			safety_x = dval;
		if (sscanf(buf, "\"safety_y\":%lf", &dval) == 1)
			safety_y = dval;

		/* read initial states */
		if (sscanf(buf, "\"vdx\":%lf", &dval) == 1)
			vdx = dval;
		if (sscanf(buf, "\"vdy\":%lf", &dval) == 1)
			vdy = dval;
		if (sscanf(buf, "\"x(0)\":%lf", &dval) == 1)
			x0 = dval;
		if (sscanf(buf, "\"y(0)\":%lf", &dval) == 1)
			y0 = dval;
		if (sscanf(buf, "\"vx(0)\":%lf", &dval) == 1)
			vx0 = dval;
		if (sscanf(buf, "\"ax(0)\":%lf", &dval) == 1)
			ax0 = dval;
		if (sscanf(buf, "\"ay(0)\":%lf", &dval) == 1)
			ay0 = dval;
		
		/* obstacles */
		if (sscanf(buf, "\"obst_x(%d,0)\":%lf", &ival, &dval) == 2)
			obsx.push_back(dval);
		if (sscanf(buf, "\"obst_y(%d,0)\":%lf", &ival, &dval) == 2)
			obsy.push_back(dval);
		if (sscanf(buf, "\"obst_vx(%d,0)\":%lf", &ival, &dval) == 2)
			obsvx.push_back(dval);
		if (sscanf(buf, "\"obst_vy(%d,0)\":%lf", &ival, &dval) == 2)
			obsvy.push_back(dval);

		memset(buf, 0, sizeof(buf));
	}

	if (!(obsx.size() == obsy.size() && obsx.size() == obsvx.size() && obsx.size() == obsvy.size())) {
		fputs("incomplete input", stderr);
		exit(1);
	}
	obs_n = obsx.size();

}

// void Data::exportSimData(class Vehicle veh, class Infrastructure inf, int it, int numsteps) {
	
// 	FILE* fd2 = NULL;
// 	if ((fd2 = fopen("../sim/data/sim.js", "w")) == NULL) {
// 		printf("\nCouldn't export data file \n");
// 	}

// 	double max_pos = 0.0;
// 	fprintf(fd2, "sim = {\n");
// 	fprintf(fd2, "x: [\n");

// 	// longitudinal position
// 	double pos_x;
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++) {
// 			pos_x = (it == -1) ? veh.init_x[k] : x[k];
// 			fprintf(fd2, "%.4f,", pos_x);
// 			if (pos_x > max_pos)
// 				max_pos = pos_x;
// 		}
// 		fputs("],\n", fd2);
// 	}

// 	for (int idVeh = 0; idVeh < obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++) {
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].x[k]);
// 			if (veh.obs[idVeh].x[k] > max_pos)
// 				max_pos = veh.obs[idVeh].x[k];
// 		}
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	// longitudinal acceleration (control)
// 	if (it != -1) {
// 		double acc_x;
// 		fputs("ax: [\n", fd2);
// 		for (int idVeh = 0; idVeh < 1; idVeh++) {
// 			fputs("\t[", fd2);
// 			for (int k = 0; k < numsteps; k++) {
// 				acc_x = (it == -1) ? veh.init_ax[k] : ux[k];
// 				fprintf(fd2, "%.4f,", acc_x);
// 			}
// 			fputs("],\n", fd2);
// 		}
// 		fputs("],\n", fd2);

// 		// lateral speed (control)
// 		double spd_y;
// 		fputs("ay: [\n", fd2);
// 		for (int idVeh = 0; idVeh < 1; idVeh++) {
// 			fputs("\t[", fd2);
// 			for (int k = 0; k < numsteps; k++) {
// 				spd_y = (it == -1) ? 0.0 : uy[k];
// 				fprintf(fd2, "%.4f,", spd_y);
// 			}
// 			fputs("],\n", fd2);
// 		}
// 		fputs("],\n", fd2);
// 	}

// 	// lateral position
// 	double pos_y;
// 	fputs("y: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++) {
// 			pos_y = (it == -1) ? veh.init_y[k] : y[k];
// 			fprintf(fd2, "%.4f,", pos_y);
// 		}
// 		fputs("],\n", fd2);
// 	}
// 	for (int idVeh = 0; idVeh < obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++)
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].y[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	// longitudinal speed
// 	double spd_x;
// 	fputs("vx: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++) {
// 			spd_x = (it == -1) ? veh.init_vx[k] : vx[k];
// 			fprintf(fd2, "%.4f,", spd_x);
// 		}
// 		fputs("],\n", fd2);
// 	}
// 	for (int idVeh = 0; idVeh < obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < numsteps; k++)
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].vx[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("id: [", fd2);
// 	for (int idVeh = 0; idVeh < 1 + obs_n; idVeh++)
// 		fprintf(fd2, "%d,", idVeh + 1);
// 	fputs("],\n", fd2);

// 	fputs("Cx: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%.4f,", 5.0);
// 	for (int idVeh = 0; idVeh < obs_n; idVeh++)
// 		fprintf(fd2, "%.4f,", 5.0);
// 	fputs("],\n", fd2);

// 	fputs("Cy: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%.4f,", 2.0);
// 	for (int idVeh = 0; idVeh < obs_n; idVeh++)
// 		fprintf(fd2, "%.4f,", 2.0);
// 	fputs("],\n", fd2);

// 	fprintf(fd2, "n:%d,\n", 1 + obs_n);
// 	fprintf(fd2, "k:%d,\n", numsteps);
// 	fprintf(fd2, "roadlength:%f,\n", max_pos);
// 	fprintf(fd2, "roadwidth:%f,\n", inf.yMax - inf.xMin);
// 	fprintf(fd2, "vdy:%f,\n", 0.0);
// 	fprintf(fd2, "Step:%f,\n", step);
// 	fprintf(fd2, "safety1:%f,\n", safety_x);
// 	fprintf(fd2, "safety2:%f,\n", safety_y);
// 	fprintf(fd2, "roadbound:%f,\n", 1.2);

// 	fprintf(fd2, "};\n\n");

// 	fclose(fd2);
// }



