#include <vector>
#include <map>

using namespace std;

#ifndef DELTASQ_H
#define DELTASQ_H

void delta_S_Q_tilde(const int pmOne, const int x, const int t, const vector<int>& start_pos,
	const double K, const double Ec_times_delta_tau, const vector<double>& G_matrix, const double lambda, const vector<int>& worm_Jx,
	const vector<int>& worm_Jt, const int Lx, const int Lt,
	const int W_x, const int W_tau, vector<double> &Qtau_tilde,
	double &delta_SA, double &delta_Qtau_tilde_val, int &debug_Q);

void delta_Qtau(const int pmOne, const int x, const int t, const int start_x, const vector<int> &worm_Jt,
	const int Lx, const int Lt, const int W_tau, int &res_accept, int &res_reject);

template<typename T>
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

#endif /* DELTASQ_H */

