#include "vehicle.h"


bool Vehicle::collision_per_veh(double here_x, double obsx, double here_y, double obsy, double here_vx, double obsvx) {
    double p = 150.0;
    double dx = here_x - obsx, dy = here_y - obsy;

    if ( exp(-sqrt( pow(dx, 4) + p*pow(dy, 4) ) / 15.0) > 0.01 )
        return true;
 
	return false;
}