#include "BAndGFunction.h"
#include <vector>
#include <math.h>
#include <numeric>
#include <iostream>

const double M_PI = 3.14159265358979323846;

using namespace std;

/*Function that fills the B_matrix with values*/
void fill_B_matrix(vector<double> &B_matrix, const int Lt) {
	for (int t = 0; t<Lt; ++t) {
		for (int t_2 = 0; t_2<Lt; ++t_2) {
			double temp = 0.0;
			for (int n = 1; n<Lt; ++n) {
				double omega = 2 * M_PI*n / (Lt);
				temp += cos(omega*(t - t_2)) / fabs(sin(omega / 2.0));
			}
			B_matrix[t_2*Lt + t] = temp;
		}
	}
	return;
}

/*Function that fills the G_matrix with values*/
void fill_G_matrix(vector<double>& G_matrix, const int Lx, const double lambda) {

	for (int i = 0; i < Lx; ++i) {
		for (int j = 0; j < Lx; ++j) {
			G_matrix[i*Lx + j] = 1.0 / Lx;
			double temp1 = 0.0;
			double temp2 = 0.0;
			for (int n = 1; n < Lx; ++n) {
				double denominator = 1.0 + 4*lambda*lambda*sin(n*M_PI /(2.0*Lx))*sin(n*M_PI /(2.0*Lx));
				temp1 += cos(n*M_PI*(i-j)/Lx)/denominator;
				temp2 += cos(n*M_PI*(i+j+1)/Lx)/denominator;
			}
			temp1 /= Lx; temp2 /= Lx;
			G_matrix[i*Lx + j] += temp1 + temp2;
		}
	}
	return;

}

