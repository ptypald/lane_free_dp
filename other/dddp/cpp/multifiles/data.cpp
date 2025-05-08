#include "data.h"

void Data::input() {

	char buf[1024];
	// const char* line;
	double dval;
	int ival;

	FILE* input;
	input = fopen("input.txt", "r");

	while (fgets(buf, sizeof(buf), input)) {
		// fprintf(stderr, "parsing: %s\n", buf);
		if (buf[0] == '\n' && buf[1] == '\n')
			break;

        if (sscanf(buf, "\"DISPLAY\":%d", &ival) == 1)
			DISPLAY = ival;
		if (sscanf(buf, "\"numsteps\":%d", &ival) == 1)
			numsteps = ival;
		if (sscanf(buf, "\"T\":%lf", &dval) == 1)
			step = dval;
		if (sscanf(buf, "\"safety_x\":%lf", &dval) == 1)
			safety_x = dval;
		if (sscanf(buf, "\"safety_y\":%lf", &dval) == 1)
			safety_y = dval;

        if (sscanf(buf, "\"UXMIN\":%lf", &dval) == 1)
			UXMIN = dval;
        if (sscanf(buf, "\"UXMAX\":%lf", &dval) == 1)
			UXMAX = dval;
        if (sscanf(buf, "\"UYMIN\":%lf", &dval) == 1)
			UYMIN = dval;
        if (sscanf(buf, "\"UYMAX\":%lf", &dval) == 1)
			UYMAX = dval;

        if (sscanf(buf, "\"VXMIN\":%lf", &dval) == 1)
			VXMIN = dval;
        if (sscanf(buf, "\"VXMAX\":%lf", &dval) == 1)
			VXMAX = dval;

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

// void Data::exportSimData(class Vehicle veh, class Infrastructure inf, int it) {
	
// 	char fileBuffer[32];
//     snprintf(fileBuffer, sizeof(char) * 32, "../sim/data/sim_%i.js", it);

// 	FILE* fd2 = NULL;
// 	if ((fd2 = fopen(fileBuffer, "w")) == NULL) {
// 		printf("\nCouldn't read data file \n");
// 	}

// 	double max_pos = 0.0;
// 	fprintf(fd2, "sim = {\n");
// 	fprintf(fd2, "x: [\n");
	
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++) {
// 			fprintf(fd2, "%.4f,", x[k]);
// 			if (x[k] > max_pos)
// 				max_pos = x[k];
// 		}
// 		fputs("],\n", fd2);
// 	}
	
// 	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++) {
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].x[k]);
// 			if (veh.obs[idVeh].x[k] > max_pos)
// 				max_pos = veh.obs[idVeh].x[k];
// 		}
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("ax: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", this->ux[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("ay: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", this->uy[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("y: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", this->y[k]);
// 		fputs("],\n", fd2);
// 	}
// 	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].y[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("vx: [\n", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", this->vx[k]);
// 		fputs("],\n", fd2);
// 	}
// 	for (int idVeh = 0; idVeh < this->obs_n; idVeh++) {
// 		fputs("\t[", fd2);
// 		for (int k = 0; k < this->numsteps; k++)
// 			fprintf(fd2, "%.4f,", veh.obs[idVeh].vx[k]);
// 		fputs("],\n", fd2);
// 	}
// 	fputs("],\n", fd2);

// 	fputs("vdx: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%.4f,", veh.getDersiredSpeedX());
// 	fputs("],\n", fd2);

// 	fputs("id: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%d,", idVeh + 1);
// 	fputs("],\n", fd2);

// 	fputs("Cx: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%.4f,", 5.0);
// 	for (int idVeh = 0; idVeh < this->obs_n; idVeh++)
// 		fprintf(fd2, "%.4f,", 5.0);
// 	fputs("],\n", fd2);

// 	fputs("Cy: [", fd2);
// 	for (int idVeh = 0; idVeh < 1; idVeh++)
// 		fprintf(fd2, "%.4f,", 2.0);
// 	for (int idVeh = 0; idVeh < this->obs_n; idVeh++)
// 		fprintf(fd2, "%.4f,", 2.0);
// 	fputs("],\n", fd2);

// 	fprintf(fd2, "n:%d,\n", 1 + this->obs_n);
// 	fprintf(fd2, "k:%d,\n", this->numsteps);
// 	fprintf(fd2, "roadlength:%f,\n", max_pos);
// 	fprintf(fd2, "roadwidth:%f,\n", inf.yMax - inf.xMin);
// 	fprintf(fd2, "vdy:%f,\n", 0.0);
// 	fprintf(fd2, "Step:%f,\n", this->step);
// 	fprintf(fd2, "safety1:%f,\n", this->safety_x);
// 	fprintf(fd2, "safety2:%f,\n", this->safety_y);
// 	fprintf(fd2, "roadbound:%f,\n", 1.2);

// 	fprintf(fd2, "};\n\n");

// 	fclose(fd2);
// }




