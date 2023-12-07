#include "LinAlg.h"
#include <iostream>


vector<double> forward_sub(vector<vector<double>> L, vector<double> b)
{
	int n = static_cast<int>(b.size());
	vector <double> y(n);

	y[0] = b[0];
	for (int i = 1; i < n; i++) {
		y[i] = b[i] - L[i][i - 1] * y[i-1];
	}

	return y;
}

vector<double> backward_sub(vector<vector<double>> U, vector<double> y)
{
	int n = static_cast<int>(y.size());
	vector <double> x(n);

	x[n-1] = y[n-1]/U[n-1][n-1];
	for (int i = n-2; i >= 0; i--) {
		x[i] = (y[i] - U[i][i + 1] * x[i + 1])/U[i][i];
	}

	return x;
}

// linear solver LU tridiag
vector<double> linear_solve_LU_no_pivoting_tridiag(vector<vector<double>> A, vector<double> b)
{
	// LU tridiagional decomp
	int n = static_cast<int>(b.size());
	vector <double> x(n);
	vector <double> y(n);
	vector <vector <double>> L(n, vector<double>(n));
	vector <vector <double>> U(n, vector<double>(n));

	for (int i = 0; i < n - 1; i++) {
		L[i][i] = 1;
		L[i + 1][i] = A[i + 1][i] / A[i][i];
		U[i][i] = A[i][i];
		U[i][i + 1] = A[i][i + 1];
		A[i+1][i+1] = A[i+1][i+1] - L[i + 1][i] * U[i][i + 1];
	}

	L[n - 1][n - 1] = 1;
	U[n - 1][n - 1] = A[n - 1][n - 1];    

	// Solve for x in Ax=b
	y = forward_sub(L, b);
	x = backward_sub(U, y);
	return x;
}


// l