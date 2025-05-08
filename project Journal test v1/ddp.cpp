#include "ddp.h"
#include <cmath>

#define Power(a,b) pow(a,b)
#define Sqrt(a) sqrt(a)

#define CONSTRAINED 1

double DDP::axUp(double axpoint) {
	return -(uxmax - axpoint);
}

double DDP::axLow(double axpoint) {
	return (uxmin - axpoint);
}

double DDP::ayUp(double aypoint) {
	return -(uymax - aypoint);
}

double DDP::ayLow(double aypoint) {
	return (uymin - aypoint);
}

double DDP::ayLB(double xPoint, double yPoint, double vyPoint, double ayPoint, double dx, double dy, double dvy, double day) {

	double T = step, offset = 1.0 + 0.001;
	double slp = 0.0, swp1 = 200.0, swp2 = 300.0, sWidth = 5.0;
	
	double result = -ayPoint - day + (2.0*Power(slp,2)*sWidth)/(exp(2.0*slp*(-swp1 + xPoint))*Power(1.0 + exp(slp*(swp1 - xPoint)),3)) - 
		(exp(slp*(swp1 - xPoint))*Power(slp,2)*sWidth)/Power(1. + exp(slp*(swp1 - xPoint)),2) - (2.0*Power(slp,2)*sWidth)/(exp(2.0*slp*(-swp2 + xPoint))*Power(1. + exp(slp*(swp2 - xPoint)),3)) + 
		(exp(slp*(swp2 - xPoint))*Power(slp,2)*sWidth)/Power(1. + exp(slp*(swp2 - xPoint)),2) - 
		(1.*dx*(-1.*exp(slp*(swp1 - xPoint)) + exp(2.0*slp*(swp1 - xPoint)))*Power(slp,3)*sWidth)/Power(1. + exp(slp*(swp1 - xPoint)),3) + 
		(dx*(-1.*exp(slp*(swp2 - xPoint)) + exp(2*slp*(swp2 - xPoint)))*Power(slp,3)*sWidth)/Power(1. + exp(slp*(swp2 - xPoint)),3) - 
		(4.*dx*exp(slp*(-swp1 + xPoint))*(-0.5 + 1.*exp(slp*(-swp1 + xPoint)))*Power(slp,3)*sWidth)/Power(1. + 1.*exp(slp*(-swp1 + xPoint)),4) + 
		(4.*dx*exp(slp*(-swp2 + xPoint))*(-0.5 + 1.*exp(slp*(-swp2 + xPoint)))*Power(slp,3)*sWidth)/Power(1. + 1.*exp(slp*(-swp2 + xPoint)),4) + 
		Klat*(-dy + (dx*exp(slp*(swp1 - xPoint))*slp*sWidth)/Power(1. + exp(slp*(swp1 - xPoint)),2) - (dx*exp(slp*(swp2 - xPoint))*slp*sWidth)/Power(1. + exp(slp*(swp2 - xPoint)),2)) + 
		(dvy - (1.*dx*(-1.*exp(slp*(swp1 - xPoint)) + exp(2*slp*(swp1 - xPoint)))*Power(slp,2)*sWidth)/Power(1. + exp(slp*(swp1 - xPoint)),3) + 
		(dx*(-1.*exp(slp*(swp2 - xPoint)) + exp(2*slp*(swp2 - xPoint)))*Power(slp,2)*sWidth)/Power(1. + exp(slp*(swp2 - xPoint)),3))*(-2.0*Sqrt(Klat) + (Klat*T)/2.) + 
		(-2.0*Sqrt(Klat) + (Klat*T)/2.)*(-((exp(slp*(swp1 - xPoint))*slp*sWidth)/Power(1. + exp(slp*(swp1 - xPoint)),2)) + 
		(exp(slp*(swp2 - xPoint))*slp*sWidth)/Power(1. + exp(slp*(swp2 - xPoint)),2) + vyPoint) + 
		Klat*(offset + sWidth/(1. + exp(slp*(swp1 - xPoint))) - sWidth/(1. + exp(slp*(swp2 - xPoint))) - yPoint);

		return result;
}

double DDP::ayUB(double yPoint, double vyPoint, double ayPoint, double dy, double dvy, double day) {

	double T = step, offset = 1.0;
	return -(-ayPoint - day - dy*Klat + dvy*(-2.0*Sqrt(Klat) + (Klat*T)/2.0) + (-2.0*Sqrt(Klat) + (Klat*T)/2.0)*vyPoint + Klat*(-yPoint + (yMax-offset)));
}

void DDP::initial_trajectory() {
    
	if (init_x.size() != numsteps + 1) {
        // fprintf(stderr, "re-adjust initial trajectories: %ld --> %d\n", init_x.size(), numsteps + 1);
		vector<double> temp_init_uy;
    
		int adjust = (int)(init_x.size() / (numsteps + 1)) + 1;
		for (int i = 0; i < init_x.size(); i++) {
			if (i % adjust == 0)
				continue;

			init_x[i] = DBL_MAX;
			init_vx[i] = DBL_MAX;
			init_y[i] = DBL_MAX;
		}

		init_x.erase(remove(init_x.begin(), init_x.end(), DBL_MAX), init_x.end());
		init_vx.erase(remove(init_vx.begin(), init_vx.end(), DBL_MAX), init_vx.end());
		init_y.erase(remove(init_y.begin(), init_y.end(), DBL_MAX), init_y.end());

		init_vy.push_back(0.0);
		// create intitial control trajectories
		for (int i = 1; i < numsteps + 1; i++) {
			init_ux.push_back((init_vx[i] - init_vx[i - 1]) / step);
			init_vy.push_back((init_y[i] - init_y[i - 1]) / step);
			init_uy.push_back((init_vy[i] - init_vy[i - 1]) / step);
		}
		init_ux.push_back(0.0);
		init_uy.push_back(0.0);
	}
}

double DDP::computeCost(vector<double> x, vector<double> y, vector<double> vx, vector<double> vy, vector<double> ax, vector<double> ay, int K) {

	double sum = 0.0;
	for (int k = 0; k < K; k++) {
		sum += 0.5*pow(ax[k], 2) + 0.5*pow(ay[k], 2) + 0.5*pow(vx[k] - vdx, 2) + + 0.5*pow(vy[k] - 0.0, 2);
	}

	return sum;
};

void DDP::initialize_parameters(Vehicle veh) {
    uxmin = veh.UXMIN; uxmax = veh.UXMAX;
    uymin = veh.UYMIN; uymax = veh.UYMAX;

    vxMin = veh.VXMIN; vxMax = veh.VXMAX;
    yMin = 0.0; yMax = 10.0;
    xMin = x0; xMax = DBL_MAX;
}

int DDP::run(Vehicle v) {
	double xPoint, yPoint, vxPoint, vyPoint, axPoint, ayPoint;

	stages *stage;
	stage = (struct Stages*)calloc(numsteps + 1, sizeof(struct Stages));

	vector<double> axInitNext, vxInitNext, xInitNext;
	vector<double> ayInitNext, vyInitNext, yInitNext;
	axInitNext.clear(); vxInitNext.clear(); xInitNext.clear(); ayInitNext.clear(); vyInitNext.clear(); yInitNext.clear();

	clock_t start, end;
    double cpu_time_used;
	start = clock();
	
	int numsteps0 = init_x.size();
	double cost = computeCost(init_x, init_y, init_vx, init_vy, init_ux, init_uy, numsteps0);
	// cout << "Initial cost: " << cost << endl;

	double w1 = 0.1, w2 = 0.1, w3 = 0.01, w4 = 0.1, w5 = 1.0;
	double T = step, vdy = 0.0;
	double dxdx, dxdy, dxdvx, dxdvy, dydx, dydy, dydvx, dydvy, dvxdx, dvxdy, dvxdvx, dvxdvy, dvydx, dvydy, dvydvx, dvydvy;
	double dxdax, dydax, dvxdax, dvydax, dxday, dyday, dvxday, dvyday;
	int it = 0;
	while (true) {
		Eigen::MatrixXf A(4,4); Eigen::MatrixXf B(4,2); Eigen::MatrixXf C(2,2);
		Eigen::MatrixXf D(2,1); Eigen::MatrixXf E(4,1);
		Eigen::MatrixXf P(4,4); Eigen::MatrixXf Q(4,1);
		Eigen::MatrixXf alpha(2,1); Eigen::MatrixXf beta(2,4);
		Eigen::MatrixXf Cstar(1,2), R(1,4);
		Eigen::MatrixXf Ieye(2,2), Ieye4(4,4);

		Ieye(0,0) = 1.0; Ieye(0,1) = 0.0; Ieye(1,0) = 0.0; Ieye(1,1) = 1.0;
		Ieye4(0,0) = 1.0; Ieye4(0,1) = 0.0; Ieye4(0,2) = 0.0; Ieye4(0,3) = 0.0;
		Ieye4(1,0) = 0.0; Ieye4(1,1) = 1.0; Ieye4(1,2) = 0.0; Ieye4(1,3) = 0.0;
		Ieye4(2,0) = 0.0; Ieye4(2,1) = 0.0; Ieye4(2,2) = 1.0; Ieye4(2,3) = 0.0;
		Ieye4(3,0) = 0.0; Ieye4(3,1) = 0.0; Ieye4(3,2) = 0.0; Ieye4(3,3) = 1.0;

		
		double mue = 0.001;
		int k = numsteps;
		while (k >= 0) {
			
			if (k == numsteps) {
				A(0,0) = 0.0; A(0,1) = 0.0; A(0,2) = 0.0; A(0,3) = 0.0; 
				A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.0; A(1,3) = 0.0;
				A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0; A(2,3) = 0.0;
				A(3,0) = 0.0; A(3,1) = 0.0; A(3,2) = 0.0; A(3,3) = 0.0;
				
				E(0,0) = 0.0; E(1,0) = 0.0; E(2,0) = 0.0; E(3,0) = 0.0;
				
				P = A; Q = E;
				P = P + (mue*Ieye4);
				
				alpha(0) = 0.0; alpha(1) = 0.0;
				beta(0,0) = 0.0; beta(0,1) = 0.0; beta(0,2) = 0.0; beta(0,3) = 0.0; 
				beta(1,0) = 0.0; beta(1,1) = 0.0; beta(1,2) = 0.0; beta(1,3) = 0.0;

				stage[k].alphaStore = alpha; stage[k].betaStore = beta;
				stage[k].PStore = P; stage[k].QStore = Q;

				k--;
				continue;
			}
			
			xPoint = init_x[k]; vxPoint = init_vx[k]; axPoint = init_ux[k];
			yPoint = init_y[k]; vyPoint = init_vy[k]; ayPoint = init_uy[k];
			
			double p11 = P(0,0), p12 = P(0,1), p13 = P(0,2), p14 = P(0,3);
			double p21 = P(1,0), p22 = P(1,1), p23 = P(1,2), p24 = P(1,3);
			double p31 = P(2,0), p32 = P(2,1), p33 = P(2,2), p34 = P(2,3);
			double p41 = P(3,0), p42 = P(3,1), p43 = P(3,2), p44 = P(3,3);

			double q1 = Q(0,0), q2 = Q(1,0), q3 = Q(2,0), q4 = Q(3,0);

			// derivatives of l
			// lxx
			dxdx = 0.0;  dxdy = 0.0; dxdvx = 0.0; dxdvy = 0.0;
			dydx = 0.0; dydy = 0.0; dydvx = 0.0; dydvy = 0.0;
			dvxdx = 0.0; dvxdy = 0.0; dvxdvx = w3; dvxdvy = 0.0;			
			dvydx = 0.0; dvydy = 0.0; dvydvx = 0.0; dvydvy = w4;

			// lx
			double xCoeff, yCoeff, vxCoeff, vyCoeff;
			xCoeff = 0.0; yCoeff = 0.0; vxCoeff = (-vdx + vxPoint)*w3; vyCoeff = (-vdy + vyPoint)*w4;

			// luu
			double daxdax, dayday, daxday, daydax;
			daxdax = w1; daxday = 0.0; daydax = 0.0; dayday = w2;

			// lu
			double axCoeff, ayCoeff;
			axCoeff = axPoint*w1; ayCoeff = ayPoint*w2;

			// lxu
			dxdax = 0.0; dydax = 0.0; dvxdax = 0.0; dvydax = 0.0;
			dxday = 0.0; dyday = 0.0; dvxday = 0.0; dvyday = 0.0;

			// derivatives of V
			// lxx + fx.Vxx.fx
			dxdx += p11;  dxdy += p12; dxdvx += p13 + p11*T; dxdvy += p14 + p12*T;
			dydx += p21; dydy += p22; dydvx += p23 + p21*T; dydvy += p24 + p22*T;
			dvxdx += p31 + p11*T; dvxdy += p32 + p12*T; dvxdvx += p33 + p13*T + T*(p31 + p11*T); dvxdvy += p34 + p14*T + T*(p32 + p12*T);
			dvydx += p41 + p21*T; dvydy += p42 + p22*T; dvydvx += p43 + p23*T + T*(p41 + p21*T); dvydvy += p44 + p24*T + T*(p42 + p22*T);
			
			// lx + fx.Vx
			xCoeff += q1; yCoeff += q2; vxCoeff += q3 + q1*T; vyCoeff += q4 + q2*T;

			// luu + fu.Vxx.fu
			daxdax += 0.5*pow(T, 2)*(p31*T + 0.5*p11*pow(T, 2)) + T*(p33*T + 0.5*p13*pow(T, 2)); 
			daxday += 0.5*pow(T, 2)*(p32*T + 0.5*p12*pow(T, 2)) + T*(p34*T + 0.5*p14*pow(T, 2)); 
			daydax += 0.5*pow(T, 2)*(p41*T + 0.5*p21*pow(T, 2)) + T*(p43*T + 0.5*p23*pow(T, 2)); 
			dayday += 0.5*pow(T, 2)*(p42*T + 0.5*p22*pow(T, 2)) + T*(p44*T + 0.5*p24*pow(T, 2));

			// lu + fu.Vx
			axCoeff += q3*T + 0.5*q1*pow(T, 2); ayCoeff += q4*T + 0.5*q2*pow(T, 2);

			// lxu + dx.Vxx.fu
			dxdax += p13*T + 0.5*p11*pow(T, 2); 
			dydax += p23*T + 0.5*p21*pow(T, 2); 
			dvxdax += 0.5*pow(T, 2)*(p31 + p11*T) + T*(p33 + p13*T); 
			dvydax = 0.5*pow(T, 2)*(p41 + p21*T) + T*(p43 + p23*T);
			
			dxday += p14*T + 0.5*p12*pow(T, 2); 
			dyday += p24*T + 0.5*p22*pow(T, 2); 
			dvxday += 0.5*pow(T, 2)*(p32 + p12*T) + T*(p34 + p14*T); 
			dvyday += 0.5*pow(T, 2)*(p42 + p22*T) + T*(p44 + p24*T);
			

			for (int o = 0; o < numObs; o++) {
				double xObs = obs[o].x[k], yObs = obs[o].y[k];
				
				dxdx += (w5*(2*Power(xObs,4) - 8.0*Power(xObs,3)*xPoint + 12.0*Power(xObs,2)*Power(xPoint,2) - 8.0*xObs*Power(xPoint,3) + 2.0*Power(xPoint,4) + 2.0*p*Power(xObs,2)*Power(yObs,2) - 4.0*p*xObs*xPoint*Power(yObs,2) + 
					2.0*p*Power(xPoint,2)*Power(yObs,2) - 4.0*p*Power(xObs,2)*yObs*yPoint + 8.0*p*xObs*xPoint*yObs*yPoint - 4.0*p*Power(xPoint,2)*yObs*yPoint + 2*p*Power(xObs,2)*Power(yPoint,2) - 
					4.0*p*xObs*xPoint*Power(yPoint,2) + 2*p*Power(xPoint,2)*Power(yPoint,2) - 2*gamma*p*Power(yObs,2)*
					Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) + 
					4.0*gamma*p*yObs*yPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2.0*gamma*p*Power(yPoint,2)*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))))/
					(2.0*exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*Power(gamma,2)*
					Power(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2),2));
				dxdy += (w5*(2*p*Power(xObs,3)*yObs - 6*p*Power(xObs,2)*xPoint*yObs + 6*p*xObs*Power(xPoint,2)*yObs - 2*p*Power(xPoint,3)*yObs + 2*Power(p,2)*xObs*Power(yObs,3) - 2*Power(p,2)*xPoint*Power(yObs,3) - 
					2.0*p*Power(xObs,3)*yPoint + 6*p*Power(xObs,2)*xPoint*yPoint - 6*p*xObs*Power(xPoint,2)*yPoint + 2*p*Power(xPoint,3)*yPoint - 6*Power(p,2)*xObs*Power(yObs,2)*yPoint + 
					6.0*Power(p,2)*xPoint*Power(yObs,2)*yPoint + 6*Power(p,2)*xObs*yObs*Power(yPoint,2) - 6*Power(p,2)*xPoint*yObs*Power(yPoint,2) - 2*Power(p,2)*xObs*Power(yPoint,3) + 
					2.0*Power(p,2)*xPoint*Power(yPoint,3) + 2*gamma*p*xObs*yObs*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2.0*gamma*p*xPoint*yObs*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2.0*gamma*p*xObs*yPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) + 
					2.0*gamma*p*xPoint*yPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))))/
					(2.0*exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*Power(gamma,2)*
					Power(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2),2));
				dydx += (w5*(2*p*Power(xObs,3)*yObs - 6*p*Power(xObs,2)*xPoint*yObs + 6*p*xObs*Power(xPoint,2)*yObs - 2*p*Power(xPoint,3)*yObs + 2*Power(p,2)*xObs*Power(yObs,3) - 2*Power(p,2)*xPoint*Power(yObs,3) - 
					2*p*Power(xObs,3)*yPoint + 6*p*Power(xObs,2)*xPoint*yPoint - 6*p*xObs*Power(xPoint,2)*yPoint + 2*p*Power(xPoint,3)*yPoint - 6*Power(p,2)*xObs*Power(yObs,2)*yPoint + 
					6*Power(p,2)*xPoint*Power(yObs,2)*yPoint + 6*Power(p,2)*xObs*yObs*Power(yPoint,2) - 6*Power(p,2)*xPoint*yObs*Power(yPoint,2) - 2*Power(p,2)*xObs*Power(yPoint,3) + 
					2*Power(p,2)*xPoint*Power(yPoint,3) + 2*gamma*p*xObs*yObs*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2*gamma*p*xPoint*yObs*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2*gamma*p*xObs*yPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) + 
					2*gamma*p*xPoint*yPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))))/
					(2.*exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*Power(gamma,2)*
					Power(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2),2));
				dydy += (w5*(2*Power(p,2)*Power(xObs,2)*Power(yObs,2) - 4*Power(p,2)*xObs*xPoint*Power(yObs,2) + 2*Power(p,2)*Power(xPoint,2)*Power(yObs,2) + 2*Power(p,3)*Power(yObs,4) - 4*Power(p,2)*Power(xObs,2)*yObs*yPoint + 
					8*Power(p,2)*xObs*xPoint*yObs*yPoint - 4*Power(p,2)*Power(xPoint,2)*yObs*yPoint - 8*Power(p,3)*Power(yObs,3)*yPoint + 2*Power(p,2)*Power(xObs,2)*Power(yPoint,2) - 
					4*Power(p,2)*xObs*xPoint*Power(yPoint,2) + 2*Power(p,2)*Power(xPoint,2)*Power(yPoint,2) + 12*Power(p,3)*Power(yObs,2)*Power(yPoint,2) - 8*Power(p,3)*yObs*Power(yPoint,3) + 
					2*Power(p,3)*Power(yPoint,4) - 2*gamma*p*Power(xObs,2)*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) + 
					4*gamma*p*xObs*xPoint*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)) - 
					2*gamma*p*Power(xPoint,2)*Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))))/
					(2.*exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*Power(gamma,2)*
					Power(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2),2));
				xCoeff += (w5*(xObs - xPoint))/(exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*gamma*
      				Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)));
				yCoeff += (w5*(p*yObs - p*yPoint))/(exp(Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2))/gamma)*gamma*
      				Sqrt(Power(xObs,2) - 2*xObs*xPoint + Power(xPoint,2) + p*Power(yObs,2) - 2*p*yObs*yPoint + p*Power(yPoint,2)));
				
				// fprintf(stderr, "xCoeff %.4f | yCoef: %.4f \n", xCoeff, yCoeff);
			}
	
			A(0,0) = dxdx; A(0,1) = dxdy; A(0,2) = dxdvx; A(0,3) = dxdvy; A(1,0) = dydx; A(1,1) = dydy; A(1,2) = dydvx; A(1,3) = dydvy; A(2,0) = dvxdx; A(2,1) = dvxdy; A(2,2) = dvxdvx; A(2,3) = dvxdvy; A(3,0) = dvydx; A(3,1) = dvydy; A(3,2) = dvydvx; A(3,3) = dvydvy;
			B(0,0) = dxdax; B(0,1) = dxday; B(1,0) = dydax; B(1,1) = dyday; B(2,0) = dvxdax; B(2,1) = dvxday; B(3,0) = dvydax; B(3,1) = dvyday;
			C(0,0) = daxdax; C(0,1) = daxday; C(1,0) = daydax; C(1,1) = dayday;
			D(0,0) = axCoeff; D(1,0) = ayCoeff;
			E(0,0) = xCoeff; E(1,0) = yCoeff; E(2,0) = vxCoeff; E(3,0) = vyCoeff;

			/* regularization */
			C = C + (mue*Ieye);
			if (C.determinant() <= 0.0) {
				// alpha(0) = fmin(-axUp(axPoint), axLow(axPoint));
				// alpha(1) = fmin(-ayUB(yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0), ayLB(xPoint, yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0, 0.0));
				// if (-ayUB(yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0) < ayLB(xPoint, yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0, 0.0)) {
				// }

				// alpha = -C.inverse()*D; beta = -C.inverse()*B.transpose();
				Eigen::MatrixXf cost_axUp;
				
				cost_axUp = (-C.inverse()*D).transpose()*C*(-C.inverse()*D) + (-C.inverse()*D).transpose()*D;
				// cout << cost_axUp << endl;
				// fprintf(stderr, "(%d, %d) C(0,0): %.4f -- C(1,1): %.4f \n", it, k, C(0,0), C(1,1));
				fprintf(stderr, "(%d, %d) C(0,1): %.8f -- C(1,0): %.8f \n", it, k, C(0,1), C(1,0));
				// cout << C.determinant() << endl;
			}
			alpha = -C.inverse()*D; beta = -C.inverse()*B.transpose();
			
			/* ~~~~~~~~~~~~~~ */

			if (1) {
				/* check if unconstrained solution violates the constraints */
				if (axUp(alpha(0) + axPoint) >= 0.0) {
					alpha(0) = -axUp(axPoint);
					beta(0,0) = 0.0; beta(0,2) = 0.0;
					beta(1,0) = 0.0; beta(1,2) = 0.0;
				}
				else if (axLow(alpha(0) + axPoint) >= 0.0) {
					alpha(0) = axLow(axPoint);
					beta(0,0) = 0.0; beta(0,2) = 0.0;
					beta(1,0) = 0.0; beta(1,2) = 0.0;
				}
				if (ayUB(yPoint, vyPoint, alpha(1) + ayPoint, 0.0, 0.0, 0.0) >= 0.0) {
					alpha(1) = -ayUB(yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0);
					beta(1,1) = 0.0; beta(1,3) = 0.0;
					beta(0,1) = 0.0; beta(1,0) = 0.0;

					Cstar(0,0) = 0.0; Cstar(0,1) = -1.0;
					R(0,0) = 0.0; R(0,1) = Klat; R(0,2) = 0.0; R(0,3) = -(-2.0*Sqrt(Klat) + (Klat*T)/2.0);
					
					beta += Cstar.transpose()*R; 
				}
				else if (ayLB(xPoint, yPoint, vyPoint, alpha(1) + ayPoint, 0.0, 0.0, 0.0, 0.0) >= 0.0) {
					alpha(1) = ayLB(xPoint, yPoint, vyPoint, ayPoint, 0.0, 0.0, 0.0, 0.0);
					beta(1,1) = 0.0; beta(1,3) = 0.0;
					beta(0,1) = 0.0; beta(1,0) = 0.0;


					Cstar(0,0) = 0.0; Cstar(0,1) = 1.0;
					R(0,0) = 0.0; R(0,1) = -Klat; R(0,2) = 0.0; R(0,3) = (-2.0*Sqrt(Klat) + (Klat*T)/2.0);
					
					beta += Cstar.transpose()*R; 
				}
			}
			/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
			
			// Vx
			Q = E + beta.transpose()*C*alpha + beta.transpose()*D + B*alpha;
			// Vxx
			P = A + beta.transpose()*C*beta + beta.transpose()*B.transpose() + B*beta;
			/* regularization */
			
			P = P + (mue*Ieye4);
			/* ~~~~~~~~~~~~~~ */
			stage[k].alphaStore = alpha; stage[k].betaStore = beta;
			stage[k].PStore = P; stage[k].QStore = Q;
			
			k--;
		}

		// Forward pass
		double xDDP, vxDDP, yDDP, vyDDP, uxDDP, uyDDP;
		double dxDDP, dyDDP, dvxDDP, dvyDDP;
		double eps = 1.0;
		Eigen::MatrixXf uDDP(2,1), uDDPrev(2,1), statesPoint(4,1);
		while(true) {
			axInitNext.clear(); vxInitNext.clear(); xInitNext.clear();
			ayInitNext.clear(); vyInitNext.clear(); yInitNext.clear();

			xDDP = init_x[0]; vxDDP = init_vx[0]; yDDP = init_y[0]; vyDDP = init_vy[0];
			xInitNext.push_back(xDDP); vxInitNext.push_back(vxDDP); yInitNext.push_back(yDDP); vyInitNext.push_back(vyDDP);

			for (int k = 0; k < numsteps; k++) {
				dxDDP = xDDP - init_x[k]; dvxDDP = vxDDP - init_vx[k];
				dyDDP = yDDP - init_y[k]; dvyDDP = vyDDP - init_vy[k];

				uDDPrev(0,0) = init_ux[k]; uDDPrev(1,0) = init_uy[k];
				statesPoint(0,0) = dxDDP; statesPoint(1,0) = dyDDP; statesPoint(2,0) = dvxDDP; statesPoint(3,0) = dvyDDP; 

				uDDP = uDDPrev + eps * (stage[k].alphaStore + stage[k].betaStore * statesPoint);
				uxDDP = uDDP(0); uyDDP = uDDP(1);

				xDDP = xDDP + vxDDP * step + 0.5 * uxDDP * pow(step, 2);
				yDDP = yDDP + vyDDP * step + 0.5 * uyDDP * pow(step, 2);
				vxDDP = vxDDP + uxDDP * step;
				vyDDP = vyDDP + uyDDP * step;

				xInitNext.push_back(xDDP); vxInitNext.push_back(vxDDP);
				yInitNext.push_back(yDDP); vyInitNext.push_back(vyDDP);
				axInitNext.push_back(uxDDP); ayInitNext.push_back(uyDDP);
			}
			break;
		}
		axInitNext.push_back(0.0); ayInitNext.push_back(0.0);

		// compute cost
		double cost = computeCost(xInitNext, yInitNext, vxInitNext, vyInitNext, axInitNext, ayInitNext, numsteps);
		// cout << "cost: " << cost << endl;

		Eigen::MatrixXf uxCurr(numsteps, 1), uxPrev(numsteps, 1);
		Eigen::MatrixXf uyCurr(numsteps, 1), uyPrev(numsteps, 1);
		for (int k = 0; k < numsteps; k++) {
			uxCurr(k, 0) = axInitNext[k];
			uxPrev(k, 0) = init_ux[k];

			uyCurr(k, 0) = ayInitNext[k];
			uyPrev(k, 0) = init_uy[k];
		}
		double norm_x = (uxCurr - uxPrev).norm(), norm_y = (uyCurr - uyPrev).norm();
		// fprintf(stderr, "The Norm is equal to: %.4f | %.4f \n", norm_x, norm_y);

		// Write trajectories for plotting
		// DDP::plotSolution(xInitNext, yInitNext, vxInitNext, vyInitNext, axInitNext, ayInitNext, it);

		// check terminal criterion, update initial trajectories
		// if (norm_x < 0.001 || it == 0) {
		if ( (norm_x < 0.1 && norm_y < 0.1) || it >= 10) {
			opt_x = xInitNext; opt_vx = vxInitNext; opt_ux = axInitNext; opt_y = yInitNext; opt_vy = vyInitNext; opt_uy = ayInitNext;
			break;
		}
		init_ux = axInitNext; init_uy = ayInitNext; init_x = xInitNext; init_y = yInitNext; init_vx = vxInitNext; init_vy = vyInitNext;

		it++;
	}

	end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(stderr, "-- DDP iterations: %d | CPU Time: %.5f \n", it, cpu_time_used);

	return 0;
}

void DDP::plotSolution(vector<double> opt_x, vector<double> opt_y, vector<double> opt_vx, vector<double> opt_vy, vector<double> opt_ux, vector<double> opt_uy, int it) {

	char fileBuffer[64];
    snprintf(fileBuffer, sizeof(char) * 64, "../outputs/ddp/journal/solutions/ddp/temp/solution_%02i.py", it);
	
	FILE *fout;
	fout = fopen(fileBuffer, "w");

	/* python file */
	fputs("solution = {\n", fout);
	
	// initial x
	fputs("'x0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_x[k]);
	fputs("],\n", fout);

	// initial y
	fputs("'y0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_y[k]);
	fputs("],\n", fout);

	// initial vx
	fputs("'vx0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_vx[k]);
	fputs("],\n", fout);

	// initial vy
	fputs("'vy0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_vy[k]);
	fputs("],\n", fout);

	// initial ux
	fputs("'ux0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_ux[k]);
	fputs("],\n", fout);

	// initial uy
	fputs("'uy0': [", fout);
	for (int k = 0; k < numsteps; k++)
		fprintf(fout, "%.4f,", init_uy[k]);
	fputs("],\n", fout);

	// vx
	fputs("'vx': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_vx[k]);
	fputs("],\n", fout);

	// vy
	fputs("'vy': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_vy[k]);
	fputs("],\n", fout);

	// x
	fputs("'x': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_x[k]);
	fputs("],\n", fout);

	// y
	fputs("'y': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_y[k]);
	fputs("],\n", fout);

	// ux
	fputs("'ux': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_ux[k]);
	fputs("],\n", fout);

	// uy
	fputs("'uy': [\n", fout);
	for (int k = 0; k < numsteps; k++) 
		fprintf(fout, "%.4f,", opt_uy[k]);
	fputs("],\n", fout);

	fprintf(fout, "'k':%d,\n", numsteps);
	fprintf(fout, "'Step':%f,\n", step);
	fprintf(fout, "'iter':%d,\n", it + 1);
	fprintf(fout, "}\n\n");


	fclose(fout);
}

