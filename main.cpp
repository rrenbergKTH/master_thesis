/**
Program that is being used in my Master's thesis that uses Monte Carlo
simulations to study phase transitions in Josephson junction chains.

@author Rasmus Renberg, October 2018
*/

#include <cstdlib>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <numeric>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <random>
#include <stdio.h>
#include <string.h>
#include "deltaSQ.h"
#include "BAndGFunction.h"
#include "deltaSxandSt.h"
#include <map>

using namespace std;

void MC_worm(const double K, const int Lx, const int Lt, const double lambda,
	const double alpha, vector<int> &Jx,
	vector<int> &Jt, vector<double> &Qtau_tilde, const vector<double> &B_matrix,
	const vector<double> &G_matrix, const int delta_tau_choice, 
	const double Ec_times_delta_tau, const vector<int> &L_directions, 
	const vector<int> &R_directions, const bool production_run);

void Ck_tilde_inv(const vector<double> &simulated_Q_tilde_values_vec, const int Lx,
	const int Lt, const double lambda, const double K, 
	const double Ec_times_delta_tau, const vector<double> &G, 
	const int N_prod_tot, double &Ck_tilde_inv);

template<typename T>
void write_to_file(const vector<T> &vector, const bool first_entry_in_row, 
	const bool start_new_file, const char* file_name);

double STD(const vector<double> &v, const double ave);
void fill_bessel_hash(const double Ec_times_beta, const double lambda, 
	const vector<double>& alpha_vec, const vector<double>& K_vec, 
	const double factor, const int remote_computer);

random_device rd;
mt19937 mt(rd());

uniform_int_distribution<int> dist_xpos;
uniform_int_distribution<int> dist_tpos;
uniform_int_distribution<int> dist_direction_1(1, 4);
uniform_int_distribution<int> dist_direction_2(0, 2);

uniform_real_distribution<double> dist_prob(0, 1);

map<double, double> bessel_list_0;
map<double, double> bessel_list_1;
map<pair<double,double>, double> bessel_list_2;
map<pair<double,double>, double> bessel_list_3;

int main(int argc, char* argv[]) {

	// System sizes and number of MC-sweeps

	double lambda = 10;
	int N_warmup = 3 * 1e4;
	int N_prod = 3 * 1e4;
	int N_block = 100;
	int remote_computer = 0; // set to 1 only if the path to bessel_list-files is other than expected (see fill_bessel_hash)

	double Ec_times_beta;
	double time_step_factor=1.0;
	int Lx = 2;
	vector<double> K_vec = {};
	vector<double> alpha_vec = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	char* output_file_name;

	// Handle input arguments
	if (argc < 5) {
		cerr << "Not enought input arguments" << endl;
		return EXIT_FAILURE;
	}
	else {
		output_file_name = argv[1];
		double K_Value;
		sscanf(argv[2], "%d", &Lx);
		double temp;
		sscanf(argv[3], "%lf", &Ec_times_beta);
	
		if (argc < 5) {
			cerr << "Not enought input arguments" << endl;
			return EXIT_FAILURE;
		}

		int k = 4;
		while (k < argc && strcmp(argv[k], "-alpha") != 0) {
			sscanf(argv[k], "%lf", &K_Value);
			K_vec.push_back(K_Value);
			k++;
		}
		if (k < argc) {
			k++;
			double alpha_value;
			alpha_vec.clear();
			while (k < argc && strcmp(argv[k], "-lambda") != 0) {
				sscanf(argv[k], "%lf", &alpha_value);
				alpha_vec.push_back(alpha_value);
				k++;
			}
			if (k < argc) {
				k++;
				sscanf(argv[k], "%lf", &lambda);
			}
			k++;
			if (k < argc) {
				k++;
				sscanf(argv[k], "%d", &N_warmup);
			}
			k++;
			if (k < argc && strcmp(argv[k], "-Nprod") == 0) {
				k++;
				sscanf(argv[k], "%d", &N_prod);
			}
			k++;
			if (k < argc && strcmp(argv[k], "-factor") == 0) {
				k++;
				sscanf(argv[k], "%lf", &time_step_factor);
			}
			k++;
			if (k < argc && strcmp(argv[k], "-remote") == 0) {
				k++;
				sscanf(argv[k], "%d", &remote_computer);
			}
		}
	}

	// read bessel files
	fill_bessel_hash(Ec_times_beta, lambda, alpha_vec, K_vec, time_step_factor, remote_computer);

	// additional parameters (not used)
	double alpha_lower = 1.49;
	double alpha_upper = 1.59;
	int extra_Nprod_factor = 20;

	// Create file and write
	write_to_file(K_vec, true, true, output_file_name);
	write_to_file(vector<double>{ Ec_times_beta }, false, false, output_file_name);
	write_to_file(vector<int>{ Lx }, false, false, output_file_name);
	write_to_file(alpha_vec, true, false, output_file_name);
	write_to_file(vector<int>{N_block}, true, false, output_file_name);
	write_to_file(vector<int>{N_warmup}, false, false, output_file_name);
	write_to_file(vector<int>{N_prod}, false, false, output_file_name);
	write_to_file(vector<int>{(int) alpha_vec.size()}, false, false, output_file_name);
	write_to_file(vector<int>{(int) K_vec.size()}, false, false, output_file_name);
	write_to_file(vector<double>{ alpha_lower, alpha_upper, 
		extra_Nprod_factor / 1.0, time_step_factor }, false, false, output_file_name);
	write_to_file(vector<double>{lambda}, false, false, output_file_name);

	// For timing of computation time
	int numberOfSweeps = K_vec.size() * alpha_vec.size()*(1 + N_block);
	int counter = 0;

	// vectors for directions on the edges
	vector<int> L_directions = { 1, 3, 4 };
	vector<int> R_directions = { 2, 3, 4 };

	// declare and set G_matrix
	vector<double> G_matrix(Lx*Lx, 0.0);
	fill_G_matrix(G_matrix, Lx, lambda);

	//  declare matrices for inverse capacitance
	vector<double> Ck_tilde_inv_vec(N_block * alpha_vec.size(), 0.0);
	vector<double> Ck_tilde_inv_mean_vec(alpha_vec.size(), 0.0);
	vector<double> Ck_tilde_inv_STD_vec(alpha_vec.size(), 0.0);

	for (int k = 0; k < K_vec.size(); ++k) { // loop over values of K
		double K = K_vec[k];

		for (int u = 0; u < alpha_vec.size(); ++u) { // loop over values of alpha
			double alpha = alpha_vec[u];

			// determine time-step and Lt
			vector<double> delta_tau_alternatives{ 1.0 / (time_step_factor*K), 
				1.0 / (time_step_factor*K*lambda), 
				2 * M_PI*alpha/time_step_factor, 
				2 * M_PI*alpha / (time_step_factor*lambda*lambda) };

			int delta_tau_choice = distance(delta_tau_alternatives.begin(), 
				min_element(delta_tau_alternatives.begin(), delta_tau_alternatives.end()));
			int Lt = ceil(Ec_times_beta / delta_tau_alternatives[delta_tau_choice]);

			cout << "choice: " << delta_tau_choice << endl;
			cout << "Lt: " << Lt << endl;
			
			double Ec_times_delta_tau = Ec_times_beta / Lt;

			// declare and fill B_matrix
			vector<double> B_matrix(Lt*Lt, 0.0);
			fill_B_matrix(B_matrix, Lt);

			// to randomize worm starting position
			decltype(dist_xpos.param()) range_x(0, Lx + 1);
			decltype(dist_tpos.param()) range_t(0, Lt - 1);
			dist_xpos.param(range_x);
			dist_tpos.param(range_t);

			// additional option (not used)
			int Nprod_tot = N_prod;

			if (alpha > alpha_lower && alpha < alpha_upper) {
				Nprod_tot *= extra_Nprod_factor;
			}

			// Define link variables
			vector<int> Jx(Lt * (Lx + 1), 0);
			vector<int> Jt(Lt * (Lx + 2), 0);
			vector<double> Qtau_tilde(Lt, 0.0);

			// Equilibration run
			for (int l = 0; l < N_warmup; ++l) {
				MC_worm(K, Lx, Lt, lambda, alpha, Jx, Jt, Qtau_tilde, B_matrix, 
					G_matrix, delta_tau_choice, Ec_times_delta_tau, L_directions, R_directions, false);
			}
			counter += 1;
			cout << ((double)counter) / numberOfSweeps * 100 << "%, Lx=" << 
				Lx << ", Ec_times_beta=" << Ec_times_beta << ", factor=" << 
				time_step_factor << ", K=" << K << ", alpha=" << alpha << 
				", lambda=" << lambda << endl;

			vector<double> Q_tilde_matrix(Nprod_tot, 0.0);
			for (int s = 0; s < N_block; ++s) { // Loop to create statistics    
													  
				for (int l = 0; l < Nprod_tot; ++l) { // Production run
							
					MC_worm(K, Lx, Lt, lambda, alpha, Jx, Jt, Qtau_tilde, B_matrix, 
						G_matrix, delta_tau_choice, Ec_times_delta_tau, L_directions, R_directions, true);
					Q_tilde_matrix[l] = accumulate(Qtau_tilde.begin(), Qtau_tilde.end(), 0.0);

				}
				counter += 1;
				cout << ((double)counter) / numberOfSweeps * 100 << "%, Lx=" << 
					Lx <<", Ec_times_beta=" << Ec_times_beta << ", factor=" << 
					time_step_factor << ", K=" << K << ", alpha=" << alpha << 
					", lambda=" << lambda << endl;

				// Compute CkTildeInv
				double CkTildeInvValue = 0.0;
						
				Ck_tilde_inv(Q_tilde_matrix, Lx, Lt, lambda, K, Ec_times_delta_tau,
					G_matrix, Nprod_tot, CkTildeInvValue);
				cout << "1/C_k=" << CkTildeInvValue << endl;
				Ck_tilde_inv_vec[s * alpha_vec.size() + u] = CkTildeInvValue;
			}

			// Get statistics (mean and standard deviation)
			vector<double> Ck_for_alpha_column = get_column(Ck_tilde_inv_vec, u, alpha_vec.size());

			Ck_tilde_inv_mean_vec[u] = accumulate(Ck_for_alpha_column.begin(),
				Ck_for_alpha_column.end(), 0.0) / Ck_for_alpha_column.size();
			Ck_tilde_inv_STD_vec[u] = STD(Ck_for_alpha_column, Ck_tilde_inv_mean_vec[u]);

		}

		// Write data to file
		cout << "Writing data to file ..." << endl;
		for (int p = 0; p < alpha_vec.size(); ++p) {
			vector<double> Ck_tilde_column = get_column(Ck_tilde_inv_vec, p, alpha_vec.size());
			write_to_file(Ck_tilde_column, true, false, output_file_name);
			
			write_to_file(vector<double> { Ck_tilde_inv_mean_vec[p] }, false, false, output_file_name);
			write_to_file(vector<double> { Ck_tilde_inv_STD_vec[p] }, false, false, output_file_name);
		}
		cout << "Done writing data to file" << endl;
	}
		
	std::cout << "DONE!";
	return EXIT_SUCCESS;
}

/*Function that generates a random closed worm in the Euclidean space-time
lattice of size Lx*Lt.*/
void MC_worm(const double K, const int Lx, const int Lt, const double lambda,
	const double alpha, vector<int> &Jx,
	vector<int> &Jt, vector<double> &Qtau_tilde, const vector<double> &B_matrix,
	const vector<double> &G_matrix, const int delta_tau_choice, 
	const double Ec_times_delta_tau, const vector<int> &L_directions, 
	const vector<int> &R_directions, const bool production_run){

	//Randomize starting position

	int startx = dist_xpos(mt);
	int startt = dist_tpos(mt);

	vector<int> start_pos = { startx, startt };
	vector<int> current_pos = { startx, startt };
	vector<int> worm_Jx(Lt * (Lx + 1), 0);
	vector<int> worm_Jt(Lt * (Lx + 2), 0);

	//Keep record of the number of steps in time and x-direction
	int W_tau = 0;
	int W_x = 0;
	int W_tau_limit = Lt;

	int debug_Q = 0;

	bool connected = false;

	while (!connected) { //While worm has not walked back to starting position

		int direction;
		if (current_pos[0] != 0 && current_pos[0] != Lx + 1) {
			direction = dist_direction_1(mt);
		}
		else if (current_pos[0] == 0) {
			direction = L_directions[dist_direction_2(mt)];
		}
		else if (current_pos[0] == Lx + 1) {
			direction = R_directions[dist_direction_2(mt)];
		}

		bool new_position = false; // False as long the worm has found no new position

		if (direction == 1) { // Try to walk in positive x-direction;

			double delta_Sx_res;
			if (current_pos[0] != 0 && current_pos[0] != Lx) { //If moving in josephson-chain
				delta_Sx_res = delta_Sx(1, current_pos[0], current_pos[1], K, 
					alpha, delta_tau_choice, lambda, Jx, Lx);
			}
			else { //else if moving to/from transmission line

				delta_Sx_res = delta_Sx_boundary(+1, current_pos[0], 
					current_pos[1], Lx, alpha, Lt, B_matrix, Jx);
			}

			double accept_rate_prefactor = 1.0;
			if (current_pos[0] == Lx) {
				accept_rate_prefactor = 4.0 / 3.0;
			}
			else if (current_pos[0] == 0) {
				accept_rate_prefactor = 3.0 / 4.0;
			}

			if (dist_prob(mt) < min(1.0, accept_rate_prefactor*exp(-delta_Sx_res))) { // If move accepted
				new_position = true;
				Jx[current_pos[1] * (Lx + 1) + current_pos[0]] += 1;
				worm_Jx[current_pos[1] * (Lx + 1) + current_pos[0]] += 1;

				current_pos[0] += 1;
				W_x += 1;

			}
		}
		else if (direction == 2) { // Try to walk in negative x-direction
			double delta_Sx_res;
			if (current_pos[0] != (Lx + 1) && current_pos[0] != 1) { // If moving in josephson-chain
				int temp_x = current_pos[0] - 1;
				delta_Sx_res = delta_Sx(-1, temp_x, current_pos[1], K, 
					alpha, delta_tau_choice, lambda, Jx, Lx);
			}
			else { //%if moving to/from transmission line
				delta_Sx_res = delta_Sx_boundary(-1, current_pos[0], 
					current_pos[1], Lx, alpha, Lt, B_matrix, Jx);
				
			}

			double accept_rate_prefactor = 1.0;
			if (current_pos[0] == 1) {
				accept_rate_prefactor = 4.0 / 3.0;
			}
			else if (current_pos[0] == Lx +1) {
				accept_rate_prefactor = 3.0 / 4.0;
			}

			if (dist_prob(mt) < min(1.0, accept_rate_prefactor*exp(-delta_Sx_res))) { // If move accepted
				new_position = true;
				current_pos[0] -= 1;

				Jx[current_pos[1] * (Lx + 1) + current_pos[0]] -= 1;
				worm_Jx[current_pos[1] * (Lx + 1) + current_pos[0]] -= 1;

				W_x -= 1;
			}

		}
		else if (direction == 3 && W_tau <= W_tau_limit) { // Try to walk in positive t-direction

			int temp_Q = 0;
			double delta_Qtau_tilde = 0.0;
			double delta_S_Q_tilde_res = 0.0;
			double delta_Stau_tot = 0.0;
			double delta_Stau_res = 0.0;
			if (current_pos[0] != 0 && current_pos[0] != Lx + 1) { // If not walking in Joesephson chain
				
				delta_S_Q_tilde(1, current_pos[0], current_pos[1], start_pos, K, 
					Ec_times_delta_tau, G_matrix, lambda, worm_Jx, worm_Jt, Lx, Lt, W_x,
					W_tau, Qtau_tilde, delta_S_Q_tilde_res, delta_Qtau_tilde, temp_Q);
				delta_Stau_res = delta_Stau(1, current_pos[0], current_pos[1], K, 
					Ec_times_delta_tau, lambda, Lx, G_matrix, Jt);
			}
			else {

				delta_S_Q_tilde(1, current_pos[0], current_pos[1], start_pos, K, 
					Ec_times_delta_tau, G_matrix, lambda, worm_Jx, worm_Jt, Lx, Lt, W_x,
					W_tau, Qtau_tilde, delta_S_Q_tilde_res, delta_Qtau_tilde, temp_Q);
				delta_Stau_res = 0.0;

			}
			delta_Stau_tot = delta_S_Q_tilde_res + delta_Stau_res;

			if (dist_prob(mt) < min(1.0, exp(-delta_Stau_tot))) { // If move accepted
				new_position = true;
				W_tau += 1;
				Jt[current_pos[1] * (Lx + 2) + current_pos[0]] += 1;
				worm_Jt[current_pos[1] * (Lx + 2) + current_pos[0]] += 1;

				Qtau_tilde[current_pos[1]] += delta_Qtau_tilde;
				current_pos[1] += 1;

				debug_Q += temp_Q;

				if (current_pos[1] == Lt) {
					current_pos[1] = 0;
				}
			}

		}
		else if (direction == 4 && W_tau >= -W_tau_limit) { // Try to walk in negative t-direction
			
			int temp_Q = 0;

			double delta_Qtau_tilde = 0.0;
			double deltaS_Q_tilde_res = 0.0;
			double delta_Stau_tot = 0.0;
			double delta_Stau_res = 0.0;
			int temp_y = current_pos[1] - 1;
			if (temp_y < 0) {
				temp_y = Lt - 1;
			}

			if (current_pos[0] != 0 && current_pos[0] != Lx + 1) { // If walking in Josephson chain

				delta_S_Q_tilde(-1, current_pos[0], temp_y, start_pos, K, 
					Ec_times_delta_tau, G_matrix, lambda, worm_Jx, worm_Jt, Lx, Lt, W_x,
					W_tau, Qtau_tilde, deltaS_Q_tilde_res, delta_Qtau_tilde, temp_Q);

				delta_Stau_res = delta_Stau(-1, current_pos[0], temp_y, K, Ec_times_delta_tau, lambda, Lx, G_matrix, Jt);
			}
			else {
				delta_S_Q_tilde(-1, current_pos[0], temp_y, start_pos, K, 
					Ec_times_delta_tau, G_matrix, lambda, worm_Jx, worm_Jt, Lx, Lt, W_x,
					W_tau, Qtau_tilde, deltaS_Q_tilde_res, delta_Qtau_tilde, temp_Q);
				delta_Stau_res = 0.0;
			}
			
			delta_Stau_tot = deltaS_Q_tilde_res + delta_Stau_res;

			if (dist_prob(mt) < min(1.0, exp(-delta_Stau_tot))) { // If move accepted

				new_position = true;
				current_pos[1] -= 1;
				if (current_pos[1] < 0) {
					current_pos[1] = Lt - 1;
				}

				Qtau_tilde[current_pos[1]] += delta_Qtau_tilde;
				debug_Q += temp_Q;
				Jt[current_pos[1] * (Lx + 2) + current_pos[0]] -= 1;
				worm_Jt[current_pos[1] * (Lx + 2) + current_pos[0]] -= 1;
				W_tau -= 1;
				
			}
		}
		if (current_pos == start_pos && new_position && W_x == 0 && W_tau == 0) {// If worm has reached starting position

			connected = true;
		}
	}

	
	/*Uncomment if you want to write the created worm to file*/
	/*if (production_run && worm_Jt != vector<int>(Lt * (Lx + 2), 0) && debug_Q!=0) {
		write_to_file(start_pos, true, true, "worm.txt");
		write_to_file(vector<int>{Lx, Lt}, false, false, "worm.txt");
		write_to_file(worm_Jx, true, false, "worm.txt");
		write_to_file(worm_Jt, true, false, "worm.txt");
		cout << debug_Q << endl;
		cout << "now";
	}*/

	return;
}

/*Function that computes the inverse capacitance*/
void Ck_tilde_inv(const vector<double> &simulated_Q_tilde_values_vec, const int Lx,
	const int Lt, const double lambda, const double K, const double Ec_times_delta_tau, const vector<double> &G_matrix,
	const int Nprod_tot, double &Ck_tilde_inv) {

	double mean = accumulate(simulated_Q_tilde_values_vec.begin(), simulated_Q_tilde_values_vec.end(), 0.0) / simulated_Q_tilde_values_vec.size();
	
	double std = STD(simulated_Q_tilde_values_vec, mean);
	double variance = std*std;

	// Compute CkTildeInv
	double factor = 1.0 / ((Lx - 1) - lambda * lambda*(G_matrix[(Lx - 1)*Lx + (Lx - 1)] + G_matrix[0] - G_matrix[(Lx - 1)*Lx] - G_matrix[(Lx - 1)]));

	Ck_tilde_inv = 1.0 - Ec_times_delta_tau / (Lt) * 1.0 / ((Lx - 1) - lambda * lambda*(G_matrix[(Lx - 1)*Lx + (Lx - 1)] + G_matrix[0] - G_matrix[(Lx - 1)*Lx] - G_matrix[(Lx - 1)]))*(variance);
	return;

}

/*Function that computes the standard deviation of a vector*/
double STD(const vector<double> &v, const double ave) {
	if (v.size() == 1) {
		return 0.0;
	}
	else {
		double res = 0.0;
		for (int i = 0; i < v.size(); ++i) {
			res += (v[i] - ave)*(v[i] - ave);
		}
		return sqrt(1.0 / (v.size() - 1) * res);
	}
}

/*Function that writes data to file_name*/
template <typename T>
void write_to_file(const vector<T> &vector, const bool first_entry_in_row,
	const bool start_new_file, const char* file_name) {
	ios_base::openmode mode;
	if (start_new_file) {
		mode = ios::trunc;
	}
	else {
		mode = ios::app;
	}
	ofstream output_file(file_name, mode);

	if (output_file.is_open()) {
		if (first_entry_in_row) {
			if (!start_new_file) {
				output_file << "\n";
			}
		}

		for (int i = 0; i < vector.size(); ++i) {
			output_file << vector[i];
			output_file << " ";
		}
		output_file.close();

	}
	else {
		cerr << "Unable to open file";
		exit(EXIT_FAILURE);
	}
	return;
}

/*Function that reads bessel-list files and puts in hash tables*/
void fill_bessel_hash(const double Ec_times_beta, const double lambda, 
	const vector<double>& alpha_vec, const vector<double>& K_vec, 
	const double factor, const int remote_computer) {
	string prefix = "";
	if (remote_computer == 1) {
		prefix = "/home/rrenberg/";
	}
	string meta_str = prefix+"bessel_list_files/bessel_list_meta_data.txt";
	string infile0_str = prefix + "bessel_list_files/bessel_list_0.txt";
	string infile1_str = prefix + "bessel_list_files/bessel_list_1.txt";
	string infile2_str = prefix + "bessel_list_files/bessel_list_2.txt";
	string infile3_str = prefix + "bessel_list_files/bessel_list_3.txt";
	string infile4_str = prefix + "bessel_list_files/bessel_list_4.txt";

	ifstream infile_meta(meta_str);
	ifstream infile0(infile0_str);
	ifstream infile1(infile1_str);
	ifstream infile2(infile2_str);
	ifstream infile3(infile3_str);
	ifstream infile4(infile4_str);

	cout << "Can read bessel files: " << infile_meta.good() << endl;
	double Ec_times_beta_read, lambda_read, alpha, K, bessel_value, factor_read;
	string line;
	getline(infile_meta, line);
	istringstream iss(line);
	vector<double> Ec_times_beta_vec_read ={istream_iterator<double>(iss), istream_iterator<double>()};

	bool Ec_times_beta_is_in_list =false;
	
	if (!(find_if(Ec_times_beta_vec_read.begin(), Ec_times_beta_vec_read.end(), [Ec_times_beta](double b) { return abs(Ec_times_beta - b) < 0.00001; }) == Ec_times_beta_vec_read.end())) {
		Ec_times_beta_is_in_list = true;
	}
	getline(infile_meta, line);
	iss.str(line);
	iss.clear();
	vector<double> factor_vec_read = { istream_iterator<double>(iss), istream_iterator<double>() };

	bool factor_is_in_list = false;

	if (!(find_if(factor_vec_read.begin(), factor_vec_read.end(), [factor](double b) { return abs(factor - b) < 0.00001; }) == factor_vec_read.end())) {
		factor_is_in_list = true;
	}
	getline(infile_meta, line);
	iss.str(line);
	iss.clear();
	vector<double> lambda_vec_read = { istream_iterator<double>(iss), istream_iterator<double>() };
	bool lambda_is_in_list = false;
	
	if (!(find_if(lambda_vec_read.begin(), lambda_vec_read.end(), [lambda](double b) { return abs(lambda - b) < 0.00001; }) == lambda_vec_read.end())) {
		lambda_is_in_list = true;
	}
	getline(infile_meta, line);
	iss.str(line);
	iss.clear();
	vector<double> alpha_vec_read = { istream_iterator<double>(iss), istream_iterator<double>() };
	bool all_alphas_are_in_list = true;
	for (double d : alpha_vec) {
		bool temp = false;
		if (!(find_if(alpha_vec_read.begin(), alpha_vec_read.end(), [d](double b) { return abs(d - b) < 0.00001; }) == alpha_vec_read.end())) {
			temp = true;
		}
		all_alphas_are_in_list = all_alphas_are_in_list && temp;
	}
	getline(infile_meta, line);
	iss.str(line);
	iss.clear();
	vector<double> K_vec_read = { istream_iterator<double>(iss), istream_iterator<double>() };
	bool all_Ks_are_in_list = true;
	for (double d : K_vec) {
		bool temp = false;
		if (!(find_if(K_vec_read.begin(), K_vec_read.end(), [d](double b) { return abs(d - b) < 0.00001; }) == K_vec_read.end())) {
			temp = true;
		}
		all_Ks_are_in_list = all_Ks_are_in_list && temp;
	}

	if (!Ec_times_beta_is_in_list || !factor_is_in_list || !lambda_is_in_list || !all_alphas_are_in_list || !all_Ks_are_in_list) {
		cout << "Error! The used parameters have not their listed bessel-values!" << endl;
		cout << "Ec_times_beta is in_list: " << Ec_times_beta_is_in_list <<
			"\ndelta_tau factor is in list: " << factor_is_in_list <<
			"\nlambda value is in list: " << lambda_is_in_list <<
			"\nall alpha values are in list: " << all_alphas_are_in_list <<
			"\nall K values are in list: " << all_Ks_are_in_list << endl;
		//exit(EXIT_FAILURE);
	}

	while (infile0 >> Ec_times_beta_read >> factor_read >> K >> bessel_value)
	{
		if (fabs(Ec_times_beta_read - Ec_times_beta) < 0.00001 && fabs(factor_read - factor) < 0.00001) {
			bessel_list_0[K] = bessel_value;
		}
	}

	while (infile1 >> Ec_times_beta_read >> factor_read >> lambda_read >> K >> bessel_value) {
		
		if (fabs(Ec_times_beta_read - Ec_times_beta) < 0.00001 && fabs(lambda_read - lambda) < 0.00001 && fabs(factor_read - factor) < 0.00001) {
			bessel_list_1[K] = bessel_value;
		}
	}
	
	while (infile2 >> Ec_times_beta_read >> factor_read >> alpha >> K >> bessel_value) {
		if (fabs(Ec_times_beta_read - Ec_times_beta) < 0.00001 && fabs(factor_read - factor) < 0.00001) {
			bessel_list_2[make_pair(alpha, K)] = bessel_value;
		}
	}
	
	while (infile3 >> Ec_times_beta_read >> factor_read >> lambda_read >> alpha >> K >> bessel_value) {
		if (fabs(Ec_times_beta_read - Ec_times_beta) < 0.00001 && fabs(lambda_read - lambda) < 0.00001 && fabs(factor_read - factor) < 0.00001) {
			bessel_list_3[make_pair(alpha,K)] = bessel_value;


		}
	}

	cout << "check:" << bessel_list_0[0.1] << ", " <<  bessel_list_1[0.1] << 
		", " << bessel_list_2[make_pair(0.1, 0.1)] << 
		", " << bessel_list_3[make_pair(0.1, 0.1)] << endl;

	return;
}
