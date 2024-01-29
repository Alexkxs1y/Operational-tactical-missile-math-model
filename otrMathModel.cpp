#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <vector>
#include "otrMathModel.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////
//						Конструктор класса						 //
///////////////////////////////////////////////////////////////////

otr::otr():J(3), stateVectorG(6), angleVector(6), d_stateVectorG(6), orientationVector(7), d_orientationVector(7), targetParams(5),
				d_targetParams(3), speedAngles(2), airVelocity(3), g(3), aeroForces(3), aeroForces_G(3), aeroTorque(3), delta(3),
				Ke(2), Kn(2), Kv(2), ssE(2), ssN(4), ssV(4), aDyn(5), bDyn(5), cDyn(2), cy(3), cz(3), mx(3), my(4), mz(4){}


///////////////////////////////////////////////////////////////////
//						Инициализация класса					 //
///////////////////////////////////////////////////////////////////

void otr::init(double &_m, vector<double> &_J, double &_d, double &_l, double &_delta_max, double &_H0, double &_V0, 
				vector<double> &_angles, vector<double> &_speedAngles, vector<double> &_angleSpeed, vector<double> &_targetPosition, vector<double> &_airVelocity, GHOST4401 &_atm,
				vector<double> &_ssParamE, vector<double> &_ssParamN, vector<double> &_ssParamV){
	m = _m; J = _J; d = _d; l = _l, delta_max = _delta_max, stateVectorG[1] = _H0; angleVector = _angles; speedAngles = _speedAngles; airVelocity = _airVelocity; atm = _atm;
	ssParamE = _ssParamE; ssParamN = _ssParamN; ssParamV = _ssParamV;

	stateVectorG[0] = 0;
	stateVectorG[2] = 0;
	
	//Такое определение начальных скоростей работает только при отсутствии ветра. 
	stateVectorG[3] = _V0 * cos(angleVector[2] - speedAngles[0]) * cos(angleVector[1] - speedAngles[1]);
	stateVectorG[4] = _V0 * sin(angleVector[2] - speedAngles[0]) * cos(angleVector[1] - speedAngles[1]);
	stateVectorG[5] = _V0 * sin(angleVector[1] - speedAngles[1]);

	for(int i = 0; i < 3; i ++){
		orientationVector[i + 4] = _angleSpeed[i];
	}
	targetParams[2] = sqrt(_targetPosition[0]*_targetPosition[0] + stateVectorG[1]*stateVectorG[1] + _targetPosition[1]*_targetPosition[1]); // Радиус вектор
	targetParams[0] = asin( - stateVectorG[1] / targetParams[2]); // Между плоскостью и линией визирования(в изначальной схеме берется больший угол)
	targetParams[1] = -atan( _targetPosition[1]/_targetPosition[0]); // Между осью Ох и линией визирования
	targetParams[3] = _targetPosition[0];
	targetParams[4] = _targetPosition[1];
	
	t = 0;
	first_init_RodGam();
	init_step_params();
}


///////////////////////////////////////////////////////////////////
//					Родриго из углов Эйлера						 //
///////////////////////////////////////////////////////////////////

void otr::first_init_RodGam(){
	orientationVector[0] = cos(angleVector[1] / 2) * cos(angleVector[0] / 2) * cos(angleVector[2] / 2) - sin(angleVector[1] / 2) * sin(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[1] = sin(angleVector[1] / 2) * cos(angleVector[0] / 2) * sin(angleVector[2] / 2) + cos(angleVector[1] / 2) * cos(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[2] = sin(angleVector[1] / 2) * cos(angleVector[2] / 2) * cos(angleVector[0] / 2) + cos(angleVector[1] / 2) * sin(angleVector[2] / 2) * sin(angleVector[0] / 2);
	orientationVector[3] = cos(angleVector[1] / 2) * sin(angleVector[2] / 2) * cos(angleVector[0] / 2) - sin(angleVector[1] / 2) * cos(angleVector[2] / 2) * sin(angleVector[0] / 2);
	double norm = sqrt(orientationVector[0]*orientationVector[0] + orientationVector[1]*orientationVector[1] +
						orientationVector[2]*orientationVector[2] + orientationVector[3]*orientationVector[3]);
	for(int i = 0; i < 4; i ++){
		orientationVector[i] /= norm;
	}
}



///////////////////////////////////////////////////////////////////
//					Дополнительные параметры					 //
///////////////////////////////////////////////////////////////////

void otr::init_extra_params(){
	Vref = sqrt( (stateVectorG[3] - airVelocity[0])*(stateVectorG[3] - airVelocity[0]) + (stateVectorG[4] - airVelocity[1])*(stateVectorG[4] - airVelocity[1]) + 
				(stateVectorG[5] - airVelocity[2])*(stateVectorG[5] - airVelocity[2])  );
	M = Vref / atm.get_a(stateVectorG[1]);
	q = Vref * Vref * atm.get_density(stateVectorG[1]) * 0.5;
	qs = q*M_PI*d*d/4;
	qsl = qs*l;
}


///////////////////////////////////////////////////////////////////
//							АДХ									 //
///////////////////////////////////////////////////////////////////

void otr::init_cx(){
	if (M < 2.16) {
		M = 2.16;
	}
	cx = 1 / (73.211 / exp(M) - 47.483 / M + 16.878);
}

void otr::init_cy(){
	double ds = 1.86 * (11.554 / exp(M) - 2.5191e-03 * M * M - 5.024 / M + 52.836e-03 * M + 4.112);
	double cy_alpha = 0.0;
	if (ds >= 0) {
		cy_alpha = sqrt(ds);
	}
	else {
		cy_alpha = 1.86 * 1.039;
	}

	double alpha_deg = speedAngles[0] * 180 / M_PI;
	double p1 = 1 / (243.84e-03 / exp(alpha_deg) + 74.309e-03);
	double p2 = log(1.9773 * alpha_deg * alpha_deg - 25.587 * alpha_deg + 83.354);
	double p3 = 18.985 * alpha_deg * alpha_deg - 375.76 * alpha_deg + 1471;
	double p4 = -51.164e-03 * alpha_deg * alpha_deg + 805.52e-03 * alpha_deg + 1.8929;
	double cy_deltav = (-p1 * 1e-06 * M * M + p2 * 1e-12 * exp(M) - p3 * 1e-06 * M - p4 * 1e-03) * 2;
	cy[1] = cy_alpha;
	cy[2] = cy_deltav;
}

void otr::init_cz(){
	cz[1]  = -cy[1]  ;
	cz[2] = -cy[2] ;
}

void otr::init_mx(){
	mx[1] = -0.005 * 0.6786 * l / Vref;

	if( abs(orientationVector[4]) > 0 ){
		mx[2] = 100  / qsl;
	} else {
		mx[2] = 0;
	}
}

void otr::init_mz(){
	double alpha_deg = speedAngles[0] * 180 / M_PI;
	mz[1] = 1.89 * (146.79e-06 * M * M - 158.98e-03 / M - 7.639e-03 * M - 68.195e-03)  * l / Vref;
	mz[2] = (-766.79e-03 / exp(M) + 438.74e-03 / M + 5.8822e-03 * M - 158.34e-03) ;
	double k1 = exp(-19.488e-03 * alpha_deg * alpha_deg - 378.62e-03 * alpha_deg + 6.7518);
	double k2 = exp(-21.234e-03 * alpha_deg * alpha_deg - 635.84e-06 * exp(alpha_deg) - 98.296e-03 * alpha_deg + 2.5938);
	mz[3] = 1.89 * sqrt(k1 * 1e-09 * M * M + k2 * 1e-06) ;
}

void otr::init_my(){
	my[1] = mz[1];
	my[2] = mz[2];
	my[3] = mz[3] ;
}

void otr::init_adch(){
	init_cx();
	init_cy();
	init_cz();
	init_mx();
	init_mz();
	init_my();
}

void otr::init_integral_adch(){
	cy[0] = cy[1] * speedAngles[0]  + cy[2] * delta[2];
	cz[0] = cz[1] * speedAngles[1] + cz[2] * delta[1];
	mx[0] = mx[1] * orientationVector[4] + mx[2] * delta[0];
	mz[0] = mz[1] * orientationVector[6]  + mz[2] * speedAngles[0] + mz[3] * delta[2];
	my[0] = my[1] * orientationVector[5] + my[2] * speedAngles[1] + my[3] * delta[1];
}


///////////////////////////////////////////////////////////////////
//					Динамические параметры						 //
///////////////////////////////////////////////////////////////////

void otr::init_aDyn(){
	double Mz_wz = mz[1] * qsl;
	double Mz_alpha = mz[2] *  qsl;
	double Mz_deltav = mz[3] *  qsl;
	double Y_alpha = cy[1] * qs;
	double Y_deltav = cy[2] * qs;
	aDyn[0] = - Mz_wz / J[2];
	aDyn[1] = - Mz_alpha / J[2];
	aDyn[2] = - Mz_deltav / J[2];
	aDyn[3] = Y_alpha / m / Vref;
	aDyn[4] = Y_deltav / m / Vref;
}

void otr::init_bDyn(){
	double My_wy = my[1] * qsl;
	double My_betta = my[2] * qsl;
	double My_deltan = my[3] * qsl;
	double Z_betta = cz[1] * qs;
	double Z_deltan = cz[2] *qs;
	bDyn[0] = - My_wy / J[1];
	bDyn[1] = - My_betta / J[1];
	bDyn[2] = - My_deltan / J[1];
	bDyn[3] = - Z_betta / m / Vref;
	bDyn[4] = - Z_deltan / m / Vref;
}

void otr::init_cDyn(){
	double Mx_wx = mx[1] * qsl;
	double M_stab = - mx[2] * qsl;
	cDyn[0] = - Mx_wx / J[0];
	cDyn[1] = - M_stab / J[0];
}

void otr::init_Dyn(){
	init_aDyn();
	init_bDyn();
	init_cDyn();
}


///////////////////////////////////////////////////////////////////
//					Динамические параметры СС 					 //
///////////////////////////////////////////////////////////////////

void otr::init_ssV(){
	ssV[0] = (aDyn[1] * aDyn[4] - aDyn[2] * aDyn[3])/ (aDyn[1] + aDyn[0] * aDyn[3]);
	ssV[1] = aDyn[2] / (aDyn[2] * aDyn[3] - aDyn[1] * aDyn[4]);
	ssV[2] = 1.0 / sqrt(aDyn[1] + aDyn[0] * aDyn[3]);
	ssV[3] = (aDyn[0] + aDyn[3]) / (2.0 * sqrt(aDyn[1] + aDyn[0] * aDyn[3]));
}

void otr::init_ssN(){
	ssN[0] = (bDyn[1] * bDyn[4] - bDyn[2] * bDyn[3])/ (bDyn[1] + bDyn[0] * bDyn[3]);
	ssN[1] = bDyn[2] / (bDyn[2] * bDyn[3] - bDyn[1] * bDyn[4]);
	ssN[2] = 1.0 / sqrt(bDyn[1] + bDyn[0] * bDyn[3]);
	ssN[3] = (bDyn[0] + bDyn[3]) / (2.0 * sqrt(bDyn[1] + bDyn[0] * bDyn[3]));
}

void otr::init_ssE(){
	ssE[0] = cDyn[1] / cDyn[0];
	ssE[1] = 1.0 / cDyn[0];
}

void otr::init_ss(){
	init_ssV();
	init_ssN();
	init_ssE();
}


///////////////////////////////////////////////////////////////////
//					Коэффициенты управления 					 //
///////////////////////////////////////////////////////////////////

void otr::init_Ke(){
	Ke[0] = (2.0 * ssE[1] * ssParamE[1] - ssParamE[0]) / (ssE[0] * ssParamE[0]);
	Ke[1] = ssE[1]/(ssE[0]*ssParamE[0]*ssParamE[0]);
}

void otr::init_Kn(){
	double D = sqrt(ssParamN[1] * ssParamN[1] * ssParamN[1] * ssParamN[1] * ssN[2] * ssN[2] -
					2.0 * ssN[3] * ssParamN[1] * ssParamN[1] * ssN[2] * ssN[1] +
					ssN[1] * ssN[1] * ssParamN[1] * ssParamN[1]);
	Kn[1] = -2.0 * ssN[2] * ( ssN[3] * ssN[1] - ssParamN[1] * ssParamN[1] * ssN[2] - D )/( ssN[0] * ssN[1] * ssN[1]);
	Kn[0] = ssParamN[0] * ( 1.0 + Kn[1] * ssN[0]) / ssN[0];
}

void otr::init_Kv(){
	double D = sqrt(ssParamV[1] * ssParamV[1] * ssParamV[1] * ssParamV[1] * ssV[2] * ssV[2] -
					2.0 * ssV[3] * ssParamV[1] * ssParamV[1] * ssV[2] * ssV[1] +
					ssV[1] * ssV[1] * ssParamV[1] * ssParamV[1]);
	Kv[1] = -2.0 * ssV[2] * ( ssV[3] * ssV[1] - ssParamV[1] * ssParamV[1] * ssV[2] - D )/( ssV[0] * ssV[1] * ssV[1]);
	Kv[0] = ssParamV[0] * ( 1.0 + Kv[1] * ssV[0]) / ssV[0];
}

void otr::init_K(){
	init_Ke();
	init_Kn();
	init_Kv();
}


///////////////////////////////////////////////////////////////////
//						Углы отклонения		 					 //
///////////////////////////////////////////////////////////////////

void otr::init_delta(){

	delta[0] = - Ke[0] * angleVector[3] - Ke[1] * angleVector[0];
	delta[1] = 10*Kn[0] * ksi_tar - Kn[1] * angleVector[4];
	delta[2] = 10*Kv[0] * phi_tar - Kv[1] * angleVector[5];

	for(int i = 0; i < 3; i++){
		if(abs(delta[i]) > delta_max){
			delta[i] = delta[i]/abs(delta[i]) * delta_max;
		}
	}
}


///////////////////////////////////////////////////////////////////
//						Силы и моменты		 					 //
///////////////////////////////////////////////////////////////////

void otr::init_g(){
	g[0] = g[2] = 0;
	g[1] = atm.get_g(stateVectorG[1]);
}

void otr::init_aeroForces(){
	aeroForces[0] = - cx * qs;
	aeroForces[1] = cy[0] * qs;
	aeroForces[2] = cz[0] * qs;
}

void otr::init_aeroForces_G(){
	double a11 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1;
	double a12 = 2 * (- orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a13 = 2 * (orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);

	double a21 = 2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a22 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1;
	double a23 = 2 * (- orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);

	double a31 = 2 * (- orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);
	double a32 = 2 * (orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);
	double a33 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[3] * orientationVector[3]) - 1;

	aeroForces_G[0] = a11 * aeroForces[0] + a12 * aeroForces[1] + a13 * aeroForces[2];
	aeroForces_G[1] = a21 * aeroForces[0] + a22 * aeroForces[1] + a23 * aeroForces[2];
	aeroForces_G[2] = a31 * aeroForces[0] + a32 * aeroForces[1] + a33 * aeroForces[2];
}

void otr::init_aeroTorque(){
	aeroTorque[0] = mx[0] * qsl;
	aeroTorque[1] = my[0] * qsl;
	aeroTorque[2] = mz[0] * qsl;
}

void otr::init_forceAndTorque(){
	init_g();
	init_aeroForces();
	init_aeroForces_G();
	init_aeroTorque();
}


///////////////////////////////////////////////////////////////////
//						Производные состояния 					 //
///////////////////////////////////////////////////////////////////

void otr::calc_d_stateVectorG(){
	for(int i = 0; i < 3; i++){
		d_stateVectorG[i] = stateVectorG[3 + i];
	}
	d_stateVectorG[3] = aeroForces_G[0] / m;
	d_stateVectorG[4] = aeroForces_G[1] / m - g[1];
	d_stateVectorG[5] = aeroForces_G[2] / m;
}

void otr::calc_d_orientationVector(){
	d_orientationVector[0] = - 0.5 * (orientationVector[4] * orientationVector[1] + orientationVector[5] * orientationVector[2] + orientationVector[6] * orientationVector[3]);
	d_orientationVector[1] = 0.5 * (orientationVector[4] * orientationVector[0] - orientationVector[5] * orientationVector[3] + orientationVector[6] * orientationVector[2]);
	d_orientationVector[2] = 0.5 * (orientationVector[4] * orientationVector[3] + orientationVector[5] * orientationVector[0] - orientationVector[6] * orientationVector[1]);
	d_orientationVector[3] = 0.5 * ( - orientationVector[4] * orientationVector[2] + orientationVector[5] * orientationVector[1] + orientationVector[6] * orientationVector[0]);

	d_orientationVector[4] = aeroTorque[0] / J[0] - (((J[2] - J[1]) * orientationVector[5] * orientationVector[6]) / J[0]);
	d_orientationVector[5] = aeroTorque[1] / J[1] - (((J[0] - J[2]) * orientationVector[4] * orientationVector[6]) / J[1]);
	d_orientationVector[6] = aeroTorque[2] / J[2] - (((J[1] - J[0]) * orientationVector[4] * orientationVector[5]) / J[2]);
}

void otr::calc_d_targetParams(){
	double V = stateVectorG[3] * cos(targetParams[1]) - stateVectorG[5] * sin(targetParams[1]);
	d_targetParams[0] = (V * sin(targetParams[0]) - stateVectorG[4] * cos(targetParams[0])) / targetParams[2];
	d_targetParams[1] = (stateVectorG[3] * sin(targetParams[1]) + stateVectorG[5] * cos(targetParams[1])) / (targetParams[2] * cos(targetParams[0]));
	//d_targetParams[2] = stateVectorG[4] * sin(targetParams[0]) - V * cos(targetParams[0]);
	targetParams[2] = sqrt( (stateVectorG[0] - targetParams[3]) * (stateVectorG[0] - targetParams[3]) + stateVectorG[1]* stateVectorG[1] +
							(stateVectorG[2] - targetParams[4]) * (stateVectorG[2] - targetParams[4]) );
	if(targetParams[2] > 500){
		phi_tar = d_targetParams[0];
		ksi_tar = d_targetParams[1];
	}
}

void otr::calc_d_dt(){
	calc_d_stateVectorG();
	calc_d_orientationVector();
	//calc_d_targetParams();
}


///////////////////////////////////////////////////////////////////
//						Интегрирование 		 					 //
///////////////////////////////////////////////////////////////////

void otr::integ_stateVectorG(double &dt){
	for(int i = 0; i < 6; i++){
		stateVectorG[i] += d_stateVectorG[i] * dt;
	}
}

void otr::integ_orientationVector(double &dt){
	for(int i = 0; i < 7; i++){
		orientationVector[i] += d_orientationVector[i] * dt;
	}
	double norm = sqrt(orientationVector[0]*orientationVector[0] + orientationVector[1]*orientationVector[1] +
						orientationVector[2]*orientationVector[2] + orientationVector[3]*orientationVector[3]);
	for(int i = 0; i < 4; i ++){
		orientationVector[i] /= norm;
	}
}

void otr::integ_targetParams(double &dt){
	for(int i = 0; i < 2; i++){
		targetParams[i] += d_targetParams[i] * dt;
	}
}

void otr::integrate(double &dt){
	integ_stateVectorG(dt);
	integ_orientationVector(dt);
	integ_targetParams(dt);
}


///////////////////////////////////////////////////////////////////
//				Определение углов и их производных 		 		 //
///////////////////////////////////////////////////////////////////

void otr::init_angleVector(){
	angleVector[0] = atan((2 * (orientationVector[0] * orientationVector[1] - orientationVector[2] * orientationVector[3])) / 
								(2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1));
	angleVector[1] = atan((2 * (orientationVector[0] * orientationVector[2] - orientationVector[3] * orientationVector[1])) / 
								(2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1));
	angleVector[2] = asin(2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]));
	
	angleVector[3] = orientationVector[4] - (tan(angleVector[2]) * (orientationVector[5] * cos(angleVector[0]) - orientationVector[6] * sin(angleVector[0])));
	angleVector[4] = ((1 / cos(angleVector[2])) * (orientationVector[5] * cos(angleVector[0]) - orientationVector[6] * sin(angleVector[0])));
	angleVector[5] = orientationVector[5] * sin(angleVector[0]) + orientationVector[6] * cos(angleVector[0]);
}

void otr::init_speedAngles(){
	double a11 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[1] * orientationVector[1]) - 1;
	double a12 = 2 * (- orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a13 = 2 * (orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);

	double a21 = 2 * (orientationVector[0] * orientationVector[3] + orientationVector[1] * orientationVector[2]);
	double a22 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[2] * orientationVector[2]) - 1;
	double a23 = 2 * (- orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);

	double a31 = 2 * (- orientationVector[0] * orientationVector[2] + orientationVector[1] * orientationVector[3]);
	double a32 = 2 * (orientationVector[0] * orientationVector[1] + orientationVector[2] * orientationVector[3]);
	double a33 = 2 * (orientationVector[0] * orientationVector[0] + orientationVector[3] * orientationVector[3]) - 1;

	double Vx = a11 * stateVectorG[3] + a21 * stateVectorG[4] + a31 * stateVectorG[5];
	double Vy = a12 * stateVectorG[3] + a22 * stateVectorG[4] + a32 * stateVectorG[5];
	double Vz = a13 * stateVectorG[3] + a23 * stateVectorG[4] + a33 * stateVectorG[5];

	speedAngles[0] = - atan(Vy/Vx);
	speedAngles[1] = asin(Vz/Vref);
}


///////////////////////////////////////////////////////////////////
//						Всё в нужном порядке	 		 		 //
///////////////////////////////////////////////////////////////////

void otr::init_step_params(){
	init_extra_params();
	init_speedAngles();
	init_angleVector();
	init_adch();
	init_Dyn();
	init_ss();
	init_K();
	calc_d_targetParams();
	init_delta();
	init_integral_adch();
	init_forceAndTorque();
	calc_d_dt();
}


///////////////////////////////////////////////////////////////////
//						Используемые функции	 		 		 //
///////////////////////////////////////////////////////////////////

void otr::update(double &dt, int & i){
	/*if(i%1000 == 0){
		cout << t << ' ' << stateVectorG[0] << ' ' << stateVectorG[1] << ' ' << stateVectorG[2] << ' ' << stateVectorG[3] << ' ' << stateVectorG[4] << ' ' << stateVectorG[5]
		<< ' ' << orientationVector[4] << ' ' << orientationVector[5] << ' ' << orientationVector[6] << ' ' << speedAngles[0] * 180.0 / M_PI << ' ' << speedAngles[1] * 180.0 / M_PI
		<< ' ' << angleVector[0] * 180.0 / M_PI << ' ' << angleVector[1] * 180.0 / M_PI << ' ' << angleVector[2] * 180.0 / M_PI
		<< ' ' << delta[0] << ' ' << delta[1] << ' ' << delta[2] << ' ' << targetParams[0]* 180.0 / M_PI << ' ' << targetParams[1]* 180.0 / M_PI << '\n';
	}*/
	
	integrate(dt);
	init_step_params();
	i++;
	t += dt;
}

vector<double> otr::get_res(double &dt, int & i){
	bool flag;
	while(abs(stateVectorG[1]) > 0.01){
		flag = true;
		if((stateVectorG[1] + dt * stateVectorG[4] < 0) && (abs(stateVectorG[1] + dt * stateVectorG[4]) > 0.01)){
			while((abs(stateVectorG[1] + dt * stateVectorG[4]) > 0.01) && (flag)){ 
				dt *= 0.5;
				if((stateVectorG[1] + dt * stateVectorG[4] > 0) && (abs(stateVectorG[1] + dt * stateVectorG[4]) > 0.01)){
					flag = false; 
				}
			}
		}
		update(dt, i);
	}

	/*cout << t << ' ' << stateVectorG[0] << ' ' << stateVectorG[1] << ' ' << stateVectorG[2] << ' ' << stateVectorG[3] << ' ' << stateVectorG[4] << ' ' << stateVectorG[5]
		<< ' ' << orientationVector[4] << ' ' << orientationVector[5] << ' ' << orientationVector[6] << ' ' << speedAngles[0] * 180.0 / M_PI << ' ' << speedAngles[1] * 180.0 / M_PI
		<< ' ' << angleVector[0] * 180.0 / M_PI << ' ' << angleVector[1] * 180.0 / M_PI << ' ' << angleVector[2] * 180.0 / M_PI
		<< ' ' << delta[0] << ' ' << delta[1] << ' ' << delta[2] << ' ' << targetParams[0]* 180.0 / M_PI << ' ' << targetParams[1]* 180.0 / M_PI << '\n';*/
		
	vector<double> res(3);
	for(int i = 0; i < 3; i++){
		res[i] = stateVectorG[i];
	}
	return res;
}