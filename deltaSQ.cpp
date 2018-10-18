#include "deltaSQ.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <numeric> 
#include <algorithm>

using namespace std;

void delta_Qtau(const int pmOne, const int x, const int t, const int start_x, const vector<int> &worm_Jt,
	const int Lx, const int Lt, const int W_tau, int &res_accept, int &res_reject) {
		
		for (int c = 1; c <= 2; ++c) { //c=1 accept scenario, c=2 reject scenario
									   // Create temporary Jx, Jt
	
			int tcoord = t;
			int q_tot=0;
	
			// Find how far the worm has wandered.
			int excess_steps_in_tau = W_tau;
			
			bool added_current_at_x = false;
			if (c == 1) { // If accept-scenario, add the time-link variable
				added_current_at_x = true;
				excess_steps_in_tau += pmOne;
				if (pmOne == 1) { // Walk to correct time coordinate
					tcoord += 1;
				}
			}
			else
				if (pmOne == -1) {
					tcoord += 1;
				}
			if (tcoord >= Lt) {
				tcoord = 0;
			}
	
			//Walk in t-direction
			int tau_direction = 0;
			if (excess_steps_in_tau > 0) {
				tau_direction = -1;
			}
			else if (excess_steps_in_tau < 0) {
				tau_direction = 1;
			}
	
			bool added_current_at_start_x = false;
			int added_current_value = 0;
			while (excess_steps_in_tau != 0) {
				if (tau_direction == -1) {
					tcoord -= 1;
					if (tcoord<0) {
						tcoord = Lt - 1;
					}
				}
	
				if (tcoord == t) {
					added_current_at_start_x = true;
					added_current_value += tau_direction;
				}
				excess_steps_in_tau += tau_direction;
				if (tau_direction == 1) {
					tcoord += 1;
					if (tcoord == Lt) {
						tcoord = 0;
					}
				}
			}
				
			//Last step, compute q on correct-time step. 
	
			int q=0;
			int i;
			i = Lx + 2;
			//Move in negative x-direction
			while (i>1) {
				q += worm_Jt[t * (Lx + 2) + i - 1] + ((added_current_at_x == true) && ((i - 1) == x))*pmOne + ((added_current_at_start_x == true) && ((i - 1) == start_x))*added_current_value;
				i -= 1;
				if (i < Lx+1 && i > 1) {
					q_tot += q;
				}
			}
	
			
			if (c == 1) {
				res_accept = q_tot;
			}
			else {
				res_reject = q_tot;
			}
		}
	return;
}

/*Function that computes deltaSA yielded from a step in the t-direction in the
Josephson chain*/
void delta_S_Q_tilde(const int pmOne, const int x, const int t, const vector<int>& start_pos,
	const double K, const double Ec_times_delta_tau, const vector<double>& G, const double lambda, const vector<int>& worm_Jx,
	const vector<int>& worm_Jt, const int Lx, const int Lt,
	const int W_x, const int W_tau, vector<double> &Qtau_tilde,
	double &delta_SA, double &delta_Qtau_tilde_val, int &debug_Q) {


	int accept_Q_res = 0;
	int reject_Q_res = 0;

	delta_Qtau(pmOne, x, t, start_pos[0], worm_Jt, Lx, Lt, W_tau, accept_Q_res, reject_Q_res);

	delta_Qtau_tilde_val = accept_Q_res - reject_Q_res;
	debug_Q = accept_Q_res - reject_Q_res;
	if (x != 0 && x != Lx + 1) {
		delta_Qtau_tilde_val -= lambda * lambda*pmOne*(G[(x - 1)*Lx + (Lx - 1)] - G[(x - 1)*Lx]);

	}

	delta_SA = Ec_times_delta_tau / (2.0)*1.0 / ((Lx - 1) - lambda * lambda*(G[(Lx - 1)*Lx + (Lx - 1)] + G[0] - G[(Lx - 1)*Lx] - G[(Lx - 1)]))*(2 * Qtau_tilde[t] * delta_Qtau_tilde_val + delta_Qtau_tilde_val * delta_Qtau_tilde_val);
	return;
}
