#include <vector>

using namespace std;

#ifndef BANDGFUNCTION_H
#define BANDGFUNCTION_H
void fill_B_matrix(vector<double>& B_matrix, const int Lt);
void fill_G_matrix(vector<double>& G_matrix, const int Lx, const double lambda);
#endif /* BANDGFUNCTION_H */

