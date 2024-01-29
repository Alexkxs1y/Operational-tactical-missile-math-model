//Работает от -2000м до 1200000м над среднем уровнем моря. Реализованы функции
//get_p(h) - давление от высоты. Па
//get_T_K(h) - термодинамическая температура от высоты. К
//get_T_m(h) - молярная температура от высоты. К
//get_T_K_from_T_m(h, Tm) - термодинамическая температура от молярной и высоты. К
//get_mol_mass(h) - молярная масса от высоты. кг/кмоль
//get_concentration(h) - концентрация молекул воздуха от высоты. 1 / м^3
//get_g(h) - ускорение свободного падения от высоты. м / с^2
//get_h_geopot(h) - геопотенциальная высота от геометрической высоты. м
//get_a(h) - скорость звука от высоты. м / с
//get_density(h) - плотность от высоты. кг / м^3
//get_mu(h) - динамическая вязкость от высоты. Применимо для высот до 90000 м. Па * с

#include <iostream>
#include <math.h>
#include "4401-81.hpp"

//Математические нужды
double polinom_for_mol_mass_low_h(double h){
	return 28.82 + 0.158 * pow(1 - 7.5 * pow(10, -8) * pow(h - 94000, 2), 0.5) - 2.479 * pow(10, -4) * pow (97000 - h, 0.5);
}

double polinom_for_mol_mass_medium_h(double h){
	if (h <= 97500){
		return polinom_for_mol_mass_low_h(97000) - 0.00012 * (h - 97000);
	} else {
		return polinom_for_mol_mass_low_h(97000) - 0.00012 * 500 - 0.0001511 * (h - 97500);
	}
}

double polinom_for_mol_mass(double h, double B0, double B1, double B2, double B3){
	return B0 + B1 * h + B2 * h * h + B3 * h * h * h;
}

double gradient_for_tempreture_T_m(double H_geopot, double H0_geopot, double T0_m, double betta_m){
	return T0_m + betta_m * (H_geopot - H0_geopot);
}

double gradient_for_temperture_T_K(double h, double h0, double T0_K, double betta){
	return T0_K + betta * (h - h0);
}

double pressure_equation_for_low_h(double p0, double H_geopot, double H0_geopot, double T_m, double T0_m, double T_K, double betta_m, double R, double g_c ){
	if (betta_m != 0){
		double lg_p = log10(p0) - (g_c * log10(T_m / T0_m) / betta_m / R);
		return pow(10, lg_p);
		
	} else{
		double lg_p = log10(p0) - 0.434294 * g_c * (H_geopot - H0_geopot) / (R * T_K);
		return pow(10, lg_p);
	}
}

double polinom_for_concentration(double h, double A0, double A1, double A2, double A3, double A4, double m){
	return (A0 + A1*h + A2*pow(h,2) + A3*pow(h,3) + A4*pow(h,4)) * pow (10, m);
}
//Конец математических нужд


void GHOST4401::init(){
	//Универсальная газовая постаянная. Дж / (К * кмоль)
	R_universal = 8314.32; 

	//Молярная масса воздуха на уровне моря. кг / кмоль
	Mol_mass = 28.964420;

	//Число Авагадро. 1 / кмоль
	N_avagadro = 602.257*pow(10,24);

	//Удельная газовая постоянная. Дж / (кг * К) 
	R = 287.05287;

	//Эмпирические коэффициенты Сатерлэнда. [S] = K, [Betta_s] = кг / (сек * м * К^0,5)
	S = 110.4;
	Betta_s = 1.458*pow(10,-6);

	//Показатель адиабаты. X_adiabat = Cp / Cv
	X_adiabat = 1.4;

	//Эффективный диаметр молекул воздуха при столкновение. м 
	Sigma_diametr = 0.365*pow(10,-9);

	//Ускорения свободного падения на среднем уровне моря. м / с^2
	g_c = 9.80665;

	//Условный радиус Земли, на котором ускорение свободного падения и вертикальный градиент ускорения на среднем уровне моря наиболее близки к истинным на широте 45°32'33". м 
	Rad_of_Earth = 6356767.0;
} 

//Начальное давление на низких высотных интервалах. Па
double GHOST4401::get_p0(int number_of_interval_h0){
	if(number_of_interval_h0 == 0)
		return 127774;

	if(number_of_interval_h0 == 1)
		return 101325;

	if(number_of_interval_h0 == 2)
		return 22632;

	if(number_of_interval_h0 == 3)
		return 5474.87;

	if(number_of_interval_h0 == 4)
		return 868.014;

	if(number_of_interval_h0 == 5)
		return 110.906;

	if(number_of_interval_h0 == 6)
		return 66.9384;

	if(number_of_interval_h0 == 7)
		return 3.95639;

	if(number_of_interval_h0 == 8)
		return 0.3634;

	if(number_of_interval_h0 == 9)
		return 0.06998;

	if(number_of_interval_h0 == 10)
		return 0.006938; 
}

//Начальная молярная температура на низких высотных интервалах. К
double GHOST4401::get_T0_m(int number_of_interval_h0){
	if(number_of_interval_h0 == 0)
		return 301.15;

	if(number_of_interval_h0 == 1)
		return 288.15;

	if(number_of_interval_h0 == 2)
		return 216.65;

	if(number_of_interval_h0 == 3)
		return 216.65;

	if(number_of_interval_h0 == 4)
		return 228.65;

	if(number_of_interval_h0 == 5)
		return 270.65;

	if(number_of_interval_h0 == 6)
		return 270.65;

	if(number_of_interval_h0 == 7)
		return 214.65;

	if(number_of_interval_h0 == 8)
		return 186.65;

	if(number_of_interval_h0 == 9)
		return 186.65;

	if(number_of_interval_h0 == 10)
		return 212.00;
} 

//Бетта_м на начальных высотных интервалах (градиент молярной температуры от геопотенциальной высота). К / м
double GHOST4401::get_betta_m(int number_of_interval_h0){
	if(number_of_interval_h0 == 0)
		return -0.0065;

	if(number_of_interval_h0 == 1)
		return -0.0065;

	if(number_of_interval_h0 == 2)
		return 0;

	if(number_of_interval_h0 == 3)
		return 0.0010;

	if(number_of_interval_h0 == 4)
		return 0.0028;

	if(number_of_interval_h0 == 5)
		return 0;

	if(number_of_interval_h0 == 6)
		return -0.0028;

	if(number_of_interval_h0 == 7)
		return -0.0020;

	if(number_of_interval_h0 == 8)
		return 0;

	if(number_of_interval_h0 == 9)
		return 0.003;

	if(number_of_interval_h0 == 10)
		return 0.011;
}
//Конец данных необходимых для апроксимации




//Начальна термодинамическая температура на больших высотных интервалах. К
double GHOST4401::get_T0_K(int number_of_interval_h0){
	if (number_of_interval_h0 == 0)
		return 334.42;

	if (number_of_interval_h0 == 1)
		return 559.6;

	if (number_of_interval_h0 == 2)
		return 695.6;

	if (number_of_interval_h0 == 3)
		return 834.4;

	if (number_of_interval_h0 == 4)
		return 941.9;

	if (number_of_interval_h0 == 5)
		return 984.65;

	if (number_of_interval_h0 == 6)
		return 995.9;

	if (number_of_interval_h0 == 7)
		return 999.9;

	if (number_of_interval_h0 == 8)
		return 1000;
}

//Бетта на больших высотных интервалах (градиент термодинамической температуры от высоты). К / м
double GHOST4401::get_betta(int number_of_interval_h0){
	if (number_of_interval_h0 == 0)
		return 0.011259;

	if (number_of_interval_h0 == 1)
		return 0.006800;

	if (number_of_interval_h0 == 2)
		return 0.003970;

	if (number_of_interval_h0 == 3)
		return 0.001750;

	if (number_of_interval_h0 == 4)
		return 0.000570;

	if (number_of_interval_h0 == 5)
		return 0.0001500;

	if (number_of_interval_h0 == 6)
		return 0.0000200;

	if (number_of_interval_h0 == 7)
		return 0.0000005;

	if (number_of_interval_h0 == 8)
		return 0;
}
//Конец данных необходимых для апроксимации


//Геопотенциальная высота в зависимости от геометрической. м
double GHOST4401::get_h_geopot(double h){
	return (Rad_of_Earth * h) / (Rad_of_Earth + h);
}

//Ускорение свободного падения от высоты. м / с^2
double GHOST4401::get_g(double h){
	return g_c * pow(Rad_of_Earth/(Rad_of_Earth + h), 2);
}

//Молярная масса от высоты. кг / кмоль
double GHOST4401::get_mol_mass (double h){
	if(h <= 94000)
		return Mol_mass;

	if(h > 94000 && h <= 97000)
		return polinom_for_mol_mass_low_h(h);

	if(h > 97000 && h <= intervals_h_for_mol_mass[0])
		return polinom_for_mol_mass_medium_h(h);

	if(h > intervals_h_for_mol_mass[0]){
		int j = 0;
		for(int i = 0; i < 6; i++){
			if(h >= intervals_h_for_mol_mass[i] && h < intervals_h_for_mol_mass[i+1])
				j = i;
		}
		return polinom_for_mol_mass(h, B[j][0], B[j][1], B[j][2], B[j][3]);
	}
}

//Термодинамическая температура(обычная) от молярной тепературы. К
double GHOST4401::get_T_K_from_T_m(double h, double T_m){
	return T_m * get_mol_mass(h) / Mol_mass;
}

//Молярная температура от высоты. К
double GHOST4401::get_T_m(double h){
	return get_T_K(h) * Mol_mass / get_mol_mass(h);
}

//Температура от высоты, К
double GHOST4401::get_T_K(double h){
	//Высота меньше 120000м и температура находится из молярной температуры, т.к. ее градиент линеен от геопотенциальной высоты
	if (h < high_h_intervals[0]){
		int j = 0;
		double h0 = 0;
		for (int i = 0; i < 11; i++){
			if ((h >= low_h_intervals[i]) && (h <= low_h_intervals[i + 1])){
				h0 = low_h_intervals[i];
				j = i;
			}
		}
		double T_m = gradient_for_tempreture_T_m(get_h_geopot(h), get_h_geopot(h0), get_T0_m(j), get_betta_m(j));
		return get_T_K_from_T_m(h, T_m);
	}
	//Высота больше 120000 и находится сразу термодинамическая температура. Т.к. в госте на таких высотах принято кусочно-линейное изменение ее от высота на данных высотах 
	if(h >= high_h_intervals[0] ){
		double h0 = 0;
		int j = 0;
		for(int i = 0; i < 9; i++){
			if(h > high_h_intervals[i] && h < high_h_intervals[i+1]){
				j = i;
				h0 = high_h_intervals[i];
			}
		}
		return gradient_for_temperture_T_K (h, h0, get_T0_K(j), get_betta(j));
	}
}

//Давление от высоты на низких высотах, полученное линеаризацией температур. Па
double GHOST4401::get_p_low_h(double h){
	double h0 = 0;
	int j = 0;
	for(int i = 0; i < 11; i++){
		if(h > low_h_intervals[i] && h < low_h_intervals[i+1]){
			h0 = low_h_intervals[i];
			j = i;
		}
	}
	return pressure_equation_for_low_h(get_p0(j), get_h_geopot(h), get_h_geopot(h0), get_T_m(h), get_T0_m(j), get_T_K(h), get_betta_m(j), R, g_c);
}

//Концентрация молекул от высоты. 1 / м^3
double GHOST4401::get_concentration(double h){
	if (h < high_h_intervals[0]){
		return N_avagadro * get_p_low_h(h) / (R_universal * get_T_K(h));
	}
	if (h > high_h_intervals[0]){
		int j = 0;
		for(int i = 0; i < 9; i++){
			if(h > intervals_h_for_concentration[i] && h < intervals_h_for_concentration[i+1])
				j = i;
		}
		return polinom_for_concentration(h, A[j][0], A[j][1], A[j][2], A[j][3], A[j][4], A[j][5]);
	}
}

//Давление от высоты. Па
double GHOST4401::get_p(double h){
	if (h < high_h_intervals[0])
		return get_p_low_h(h);
		
	if (h >= high_h_intervals[0])
		return get_concentration(h) * R_universal * get_T_K(h) / N_avagadro;
} 

//Скорость звука от высоты. м/с
double GHOST4401::get_a(double h){
	return std::sqrt(X_adiabat * R_universal * get_T_K(h) / get_mol_mass(h));
}

//Плотность от высоты. кг / м^3
double GHOST4401::get_density(double h){
	return  get_p(h) * get_mol_mass(h) / (R_universal * get_T_K(h));
}

//Динамическая вязкость от высоты. Применимо для высот до 90000 м. Па * с
double GHOST4401::get_mu(double h){
	return Betta_s * std::sqrt(get_T_K(h)) / (get_T_K(h) + S);
}


