#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <tuple>
using namespace std;

// forward substitution for lower triangular matrix
vector<double> forward_sub(vector<vector<double>> L, vector<double> b);

// backward substitution for upper triangular matrix
vector<double> backward_sub(vector<vector<double>> U, vector<double> b);

// linear solver LU tridiag
vector<double> linear_solve_LU_no_pivoting_tridiag(vector<vector<double>> A, vector<double> b);

#endif // !LINALG_H
