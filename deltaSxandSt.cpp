
#include "deltaSxandSt.h"
#include <iostream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

/*Function computes deltaS_t yielded by a step in the t-direction in the
Josephson-chain*/
double delta_Stau(const int pmOne, const int x, const int t, const double K, const double Ec_times_delta_tau, const double lambda, const int Lx,
	const vector<double> &G_matrix, const vector<int> &Jt) {
	double Gxx = G_matrix[(x - 1)*Lx + (x - 1)];
	double Gres = 0.0;
	for (int i = 1; i<Lx + 1; ++i) {
		if (i != x) {

			Gres += G_matrix[(x - 1)*Lx + i - 1] * Jt[t*(Lx + 2) + i];
		}
	}
	
	return (Ec_times_delta_tau*lambda*lambda/2.0)*Gxx*(1.0 + pmOne * 2.0*Jt[t*(Lx + 2) + x]) + 2.0*pmOne*(Ec_times_delta_tau*lambda*lambda / 2.0)*Gres;
}
/*Function computes deltaS_x yielded by a step in the x-direction in the
Josephson-chain*/
double delta_Sx(const int pmOne, const int x, const int t, const double K, const double alpha, const int delta_tau_choice, const double lambda, const vector<int> &Jx, const int Lx) {
	
	double bessel_value=0.0;
	
	//cout << delta_tau_choice;
	switch (delta_tau_choice)
	{
	case 0: bessel_value = bessel_list_0[K];
		break;
	case 1: bessel_value = bessel_list_1[K];
		break;
	case 2: bessel_value = bessel_list_2[make_pair(alpha,K)];
		break;
	case 3: bessel_value = bessel_list_3[make_pair(alpha,K)];
		break;
	}

	return -bessel_value*(pmOne*2.0*Jx[t*(Lx + 1) + x] + 1.0);
}

/*Function computes deltaS_x yielded by a step in the x-direction in to or out of
the Josephson-chain*/
double delta_Sx_boundary(const int pmOne, const int x, const int t, const int Lx, const double alpha, const int Lt,
	const vector<double> &B_matrix, const vector<int> &Jx) {
	
	int column;
	/* Find the correct Jx-links, those on either the left or right side of the
	Josephson chain*/
	if (pmOne == 1) {
		if (x == 0) {
			column = 0;
		}
		else {
			column = Lx;
		}
	}
	else {
		if (x == 1) {
			column = 0;
		}
		else {
			column = Lx;
		}
	}

	// Now, compute deltaSx
	double preConst = 2.0*M_PI*alpha / Lt;
	double temp = 0.0;
	
	for (int s = 0; s<Lt; ++s) {
		
		temp += pmOne * Jx[s*(Lx + 1) + column] * (B_matrix[s*Lt + t]);
	}
	temp += B_matrix[0]/2.0;
	

	return preConst * temp;
}

