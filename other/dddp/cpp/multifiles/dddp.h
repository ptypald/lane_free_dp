#ifndef DDDP_H
#define DDDP_H

#include "controller.h"
#include <vector>

class DDDP : public Controller {
    public:
        vector<double> init_vx, init_x, init_y;
        // corridor parameters
        double dux0, duy0, cvx0, cxm, cy0;
        double cx, cy, cvx;

        DDDP(int K, double T, string m) {
		numsteps = K;
		step = T;
		method = m;
	}

    void feasible_domain(int it);
    void initial_path();
    void execute(Vehicle v);

};

#endif