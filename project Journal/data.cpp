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
		if (sscanf(buf, "\"timeHorizon\":%lf", &dval) == 1)
			this->timeHorizon = dval;
		if (sscanf(buf, "\"safety_x\":%lf", &dval) == 1)
			safety_x = dval;
		if (sscanf(buf, "\"safety_y\":%lf", &dval) == 1)
			safety_y = dval;
		if (sscanf(buf, "\"DISPLAY\":%d", &ival) == 1)
			DISPLAY = ival;

		// road parameters
		if (sscanf(buf, "\"yMin\":%lf", &dval) == 1)
			yMin = dval;
		if (sscanf(buf, "\"yMax\":%lf", &dval) == 1)
			yMax = dval;
		if (sscanf(buf, "\"xMin\":%lf", &dval) == 1)
			xMin = dval;
		if (sscanf(buf, "\"xMax\":%lf", &dval) == 1)
			xMax = dval;

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


