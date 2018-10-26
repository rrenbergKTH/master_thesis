#include <vector>
#include <map>
//#include "getColumn.tpp"

using namespace std;

#ifndef DELTASXANDST_H
#define DELTASXANDST_H

extern map<double, double> bessel_list_0;
extern map<double, double> bessel_list_1;
extern map<pair<double,double>, double> bessel_list_2;
extern map<pair<double,double>, double> bessel_list_3;

double delta_Stau(const int pmOne, const int x, const int t, const double K, const double Ec_times_delta_tau, const double lambda, const int Lx,
	const vector<double> &G_matrix, const vector<int> &Jt);
double delta_Sx(const int pmOne, const int x, const int t, const double K, const double alpha, const int delta_tau_choice, const double lambda, const vector<int> &Jx, const int Lx);
double delta_Sx_boundary(const int pmOne, const int x, const int t, const int Lx, const double alpha, const int Lt,
	const vector<double> &B_matrix, const vector<int> &Jx);

/*Function that returns the column of a matrix of dimension xlength x ylength*/
template<typename T>
vector<T> get_column(const vector<T> &matrix, const int column, const int x_length) {
	int y_length = matrix.size() / x_length;
	vector<T> column_vec(y_length, 0);
	for (int i = 0; i<y_length; ++i) {
		column_vec[i] = matrix[i*x_length + column];
	}
	return column_vec;
}


#endif /* DELTASXANDST_H */

#pragma once
