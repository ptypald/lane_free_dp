#include "controllers.h"


void Controllers::obst_prediction(vector<double> x0, vector<double> y0, vector<double> vx0){

    for (int i = 0; i < numObs; i++) {	
        obs[i].x.push_back(x0[i]); obs[i].y.push_back(y0[i]);
        obs[i].vx.push_back(vx0[i]);
    }

    // allocate obs struct 
    for (int i = 0; i < numObs; i++) {
        for (int k = 0; k < numsteps; k++) {
            obs[i].x.push_back(obs[i].x[k] + obs[i].vx[k]*step);
            obs[i].y.push_back(obs[i].y[k]);
            obs[i].vx.push_back(obs[i].vx[k]);
        }
    }
    fprintf(stderr, "-- Obstacles path created!\n");

}

void Controllers::exportSimData(int it) {
	
	FILE* fd2 = NULL;
	if ((fd2 = fopen("../sim/data/sim.js", "w")) == NULL) {
		printf("\nCouldn't export data file \n");
	}

	double max_pos = 0.0;
	fprintf(fd2, "sim = {\n");
	fprintf(fd2, "x: [\n");

	// longitudinal position
	double pos_x;
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			pos_x = (it == -1) ? init_x[k] : opt_x[k];
			fprintf(fd2, "%.4f,", pos_x);
			if (pos_x > max_pos)
				max_pos = pos_x;
		}
		fputs("],\n", fd2);
	}

	for (int idVeh = 0; idVeh < numObs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			fprintf(fd2, "%.4f,", obs[idVeh].x[k]);
			if (obs[idVeh].x[k] > max_pos)
				max_pos = obs[idVeh].x[k];
		}
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	// lateral position
	double pos_y;
	fputs("y: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			pos_y = (it == -1) ? init_y[k] : opt_y[k];
			fprintf(fd2, "%.4f,", pos_y);
		}
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < numObs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++)
			fprintf(fd2, "%.4f,", obs[idVeh].y[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	// longitudinal speed
	double spd_x;
	fputs("vx: [\n", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++) {
			spd_x = (it == -1) ? init_vx[k] : opt_vx[k];
			fprintf(fd2, "%.4f,", spd_x);
		}
		fputs("],\n", fd2);
	}
	for (int idVeh = 0; idVeh < numObs; idVeh++) {
		fputs("\t[", fd2);
		for (int k = 0; k < numsteps; k++)
			fprintf(fd2, "%.4f,", obs[idVeh].vx[k]);
		fputs("],\n", fd2);
	}
	fputs("],\n", fd2);

	fputs("id: [", fd2);
	for (int idVeh = 0; idVeh < 1 + numObs; idVeh++)
		fprintf(fd2, "%d,", idVeh + 1);
	fputs("],\n", fd2);

	fputs("Cx: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	for (int idVeh = 0; idVeh < numObs; idVeh++)
		fprintf(fd2, "%.4f,", 5.0);
	fputs("],\n", fd2);

	fputs("Cy: [", fd2);
	for (int idVeh = 0; idVeh < 1; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	for (int idVeh = 0; idVeh < numObs; idVeh++)
		fprintf(fd2, "%.4f,", 2.0);
	fputs("],\n", fd2);

	fprintf(fd2, "n:%d,\n", 1 + numObs);
	fprintf(fd2, "k:%d,\n", numsteps);
	fprintf(fd2, "roadlength:%f,\n", max_pos);
	fprintf(fd2, "roadwidth:%f,\n", 10.0);
	fprintf(fd2, "vdy:%f,\n", 0.0);
	fprintf(fd2, "Step:%f,\n", step);
	fprintf(fd2, "roadbound:%f,\n", 1.2);
	fprintf(fd2, "};\n\n");

	fclose(fd2);
}


