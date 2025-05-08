/* 3 state variables -> x, y, vx */
/* controls -> long. acceleration and lateral speed */

#include <stdio.h>
#include <time.h>
#include <iostream>

#include "data.h"
#include "infrastructure.h"
#include "vehicle.h"
#include "controller.h"
#include "dddp.h"

using namespace std;
 
int main() {

	clock_t tic, toc;
	double cpu_time;
	
	tic = clock();
	/* start procedure */

	Data data;
	data.input();
	
	// Infrastructure road;

	Vehicle veh(data.x0, data.y0, data.vx0);
	veh.setDersiredSpeedX(data.vdx);
	
	veh.obs_n = data.obs_n;
	for (int i = 0; i < veh.obs_n; i++) {
		veh.obsx0.push_back(data.obsx[i]); veh.obsy0.push_back(data.obsy[i]);
		veh.obsvx0.push_back(data.obsvx[i]); veh.obsvy0.push_back(data.obsvy[i]);
	}

	veh.setUXmin(data.UXMIN); veh.setUXmax(data.UXMAX);
	veh.setUYmin(data.UYMIN); veh.setUYmax(data.UYMAX);
	veh.setVXmin(data.VXMIN); veh.setVXmax(data.VXMAX);
	veh.step = data.step; veh.numsteps = data.numsteps;

  
	DDDP dddp(data.numsteps, data.step, "dddp");
	dddp.dux0 = 1.0; dddp.duy0 = 1.0; dddp.cvx0 = 3.0; dddp.cxm = 5.0; dddp.cy0 = 10.0;
	dddp.max_iter = 0;
	dddp.x0 = veh.getX0(); dddp.y0 = veh.getY0(); dddp.vx0 = veh.getVX0(); 
	dddp.vdx = veh.getDersiredSpeedX(); dddp.vdy = veh.getDersiredSpeedY();
	dddp.uxLB = veh.getUXmin(); dddp.uxUB = veh.getUXmax();
	dddp.uyLB = veh.getUYmin(); dddp.uyUB = veh.getUYmax();
	dddp.vxLB = veh.getVXmin(); dddp.vxUB = veh.getVXmax();
	// 
	dddp.roadxMax = 10000.0; dddp.roadxMin = 0.0; dddp.roadyMax = 10.0; dddp.roadyMin = 0.0;

	dddp.initial_path();
	dddp.execute(veh);


	/* end procedure */
	toc = clock();
	cpu_time = (double)(toc - tic) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- Total CPU Time used: %.6f --\n", cpu_time);

	return 0;
}
