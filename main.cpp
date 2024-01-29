#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include "otrMathModel.hpp"

int main(){
	double m = 1500;
	vector<double> J = {170, 640, 640};
	double d = 0.98;
	double l = 7;
	double delta_max = 20 * M_PI / 180;
	double H0 = 9550;
	double V0 = 1210;
	vector<double> angles = {0, 0, -60 * M_PI / 180};
	vector<double> speedAngles = {0, 0, 0};
	vector<double> angleSpeed = {0.01, 0.01, 0.01};
	vector<double> targetPosition = {6822.66, 0};
	vector<double> airVelocity = {0, 0, 0};
	GHOST4401 atm;
	atm.init();
	
	vector<double> ssParamE = {0.01, 0.35};
	vector<double> ssParamN = {0.95, 0.35};
	vector<double> ssParamV = {0.95, 0.35};
	int iter = 0;

	otr la;
	la.init(m, J, d, l, delta_max, H0, V0, angles, speedAngles, angleSpeed, targetPosition, airVelocity, atm, ssParamE, ssParamN, ssParamV);
	double dt = 0.0001;
	//la.get_res(dt, iter);
	bool flag = true;
	double step = 200;
	targetPosition[0] = 5000;
	targetPosition[1] = 0;
	double miss = 0;
	double counter;
	for(int j = 0; j < 37; j++){
		flag = true;
		step = 100;
		targetPosition[0] = 5000;
		targetPosition[1] = 0;
		counter = 1;
		miss = 0;
		while(flag){

			la.init(m, J, d, l, delta_max, H0, V0, angles, speedAngles, angleSpeed, targetPosition, airVelocity, atm, ssParamE, ssParamN, ssParamV);
			vector<double> res = la.get_res(dt);
			miss =  sqrt( (res[0] - targetPosition[0]) * (res[0] - targetPosition[0]) +  (res[2] - targetPosition[1]) * (res[2] - targetPosition[1]) );
			if(abs( miss - 15 ) < 0.5){
				flag = false;
				cout << targetPosition[0] << ' ' << targetPosition[1] << ' ' << res[0] << ' ' << res[2] << '\n';
			} else {
				if(miss < 15){
					targetPosition[0] += step * cos(M_PI * double(j) / 36);
					targetPosition[1] += step * sin(M_PI * double(j) / 36);
				} else {
					targetPosition[0] -= step * cos(M_PI * double(j) / 36);
					targetPosition[1] -= step * sin(M_PI * double(j) / 36);
					step *= 0.5;
				}
			}
		}
	}

	return 1;
}