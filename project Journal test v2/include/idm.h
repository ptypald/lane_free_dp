#ifndef IDM_H
#define IDM_H

#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <float.h>
#include <algorithm>

#include "controllers.h"

using namespace std;

class IDM : public Controllers {
    
    public:
        double aMax = 1.5, b = 2.0, d = 4, s0 = 5.0;
        vector<double> x, y, vx, ax;

        double xl, vl, xf, vf;
        int hasLeader, hasFollower, leaderID, followerID;

    void idm_run();
    double idm_accel(double x, double xl, double v, double vl, double v0, double ts, int leaderID);
};

#endif