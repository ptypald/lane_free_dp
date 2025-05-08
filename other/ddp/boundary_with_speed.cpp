// This code calculates the boundary's value and is based from Karteek's code. This implementation is a simplified form of an internal function within the TrafficFluid-Sim application.
double boundary_value(std::vector<double>& lim, std::vector<double>& slope, std::vector<double>& offset, double pos_x, double veh_speed_x, double* boundary_speed) {	
	// based on Sec. 3.2 on IBC paper:
	// lim: the consecutive lateral levels: y_1, y_2, y_3, ...
	// slope: the consecutive slope values: s_1, s_2, ...
	// offset: the consecutive offset values: o_1, o_2, ...
	double boundary_val_tmp = 0;
	double boundary_speed_tmp = 0;
	
	double mid_height = 0.5; // the n_p parameter relevant for assymetric boundary changes, as discussed at the end of Sec. 3.1
	double direction = 1; // if direction = -1, this would work for opposite direction paths (e.g., as in IBC for direction b)
	int size_points = lim.size();
	for (size_t i = 0; i < size_points - 1; i++) {
		if (pos_x * direction < offset[i] * direction) {
			boundary_val_tmp += mid_height * (lim[i + 1] - lim[i]) * tanh(slope[i] * direction * (pos_x - offset[i]));
			if (boundary_speed != nullptr) {
				boundary_speed_tmp += mid_height * (lim[i + 1] - lim[i]) * pow(cosh(slope[i] * direction * (long_pos - offset[i])), -2) * slope[i] * veh_speed_x;
			}
		}
		else {
			boundary_val_tmp += (1 - mid_height) * (lim[i + 1] - lim[i]) * tanh((mid_height / (1 - mid_height)) * slope[i] * direction * (pos_x - offset[i]));					
			if (boundary_speed != nullptr) {
				boundary_speed_tmp += mid_height * slope[i] * veh_speed_x * (lim[i + 1] - lim[i]) * pow(cosh((mid_height / (1 - mid_height)) * slope[i] * direction * (long_pos - offset[i])), -2);
			}
		}
	}


	if (boundary_speed != nullptr) {
		*boundary_speed = boundary_speed_tmp;
	}
	
	return boundary_val_tmp + mid_height * (lim[size_points - 1] + lim[0]);
	
}
