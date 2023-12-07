#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <numeric>

#include "BlackScholes.hpp"
#include "FiniteDifference.hpp"

using namespace std;

void FiniteDifferenceDaoOptionPricer::discretize(double alpha_tmp) {
	// Domain parameter
	x_compute = log(option_base.S / option_base.K);
	x_left = log(option_base.B / option_base.K);
	tau_final = option_base.T * option_base.v * option_base.v / 2.0;
	delta_tau = tau_final / M;
	N_left = int((x_compute - x_left) / sqrt(delta_tau / alpha_tmp));
	delta_x = (x_compute - x_left) / N_left;
	x_right = log(option_base.S / option_base.K)
			+ (option_base.r - option_base.q - 0.5 * pow(option_base.v, 2))
			* option_base.T + 3 * option_base.v * sqrt(option_base.T);
	N_right = ceil((x_right - x_compute) / delta_x);
	x_right = x_compute + N_right * delta_x;
	N = N_left + N_right;
	alpha = delta_tau / delta_x / delta_x;

	return;
}

// ---------------------------------------------------
// forward Euler
// ---------------------------------------------------
void FiniteDifferenceDaoOptionPricer::computeDomainFE() {
	domain.resize(M+1, std::vector<double>(N+1));
	for (int i = 0; i < N + 1; ++i) {
		double x_i = x_left + i * delta_x;
		domain[0][i] = option_base.K * exp(a * x_i) * max(exp(x_i) - 1, 0.0);
	}
	for (int i = 0; i < M + 1; ++i) {
		double tau_i = i * delta_tau;
		double tmp = tau_i / option_base.v / option_base.v;
		domain[i][0] = 0;
		domain[i][N] = option_base.K * exp(a * x_right + b * tau_i)
				 	* (exp(x_right - 2 * option_base.q * tmp)
				 	- exp(-2 * option_base.r * tmp));
	}

	//calculte for m = 1:M
	for (int m = 1; m <= M; ++m) {
		for (int j = 1; j <= N - 1; ++j) {
			domain[m][j] = alpha * domain[m - 1][j + 1]
						+ (1 - 2 * alpha)*domain[m - 1][j]
						+ alpha*domain[m - 1][j - 1];
		}
	}
	return;
}

// ---------------------------------------------------
// backward Euler with LU
// ---------------------------------------------------
void FiniteDifferenceDaoOptionPricer::initDomain() {
	domain.resize(M+1, std::vector<double>(N+1));
	A.resize(N-1, vector<double>(N-1));

	// set up domain
	for (int n = 0; n < N; ++n) {
		domain[0][n] = f(x_left + n * delta_x);
	}
	for (int m = 0; m < M+1; ++m) {
		domain[m][N] = g_right(m * delta_tau);
		domain[m][0] = 0;
	}
}

void FiniteDifferenceDaoOptionPricer::computeDomainBE() {
	initDomain();

	// set up A
	for (int i = 1; i < N - 2; ++i) {
		A[i][i - 1] = -alpha;
		A[i][i] = 1 + 2 * alpha;
		A[i][i + 1] = -alpha;
	}
	A[0][1] = -alpha;
	A[N-2][N-3] = -alpha;

	A[0][0] = 1 + 2 * alpha;
	A[N-2][N-2] = 1 + 2 * alpha;

	// update domain
	for (int m = 1; m < M + 1; ++m) {
		vector<double> b_vec(N-1);
		for (int i = 1; i < N - 2; ++i) {
			b_vec[i] = domain[m - 1][i + 1];
		}
		b_vec[0] = domain[m - 1][1] + alpha * domain[m][0];
		b_vec[N - 2] = domain[m - 1][N - 1] + alpha * domain[m][N];

        // Solve for Ax=b
		auto x = linear_solve_LU_no_pivoting_tridiag(A, b_vec);
		for (int n = 0; n < N-1; ++n) {
			domain[m][n+1] = x[n];
		}
	}

	return;
}


// ---------------------------------------------------
// Crank Nicolson
// ---------------------------------------------------
void FiniteDifferenceDaoOptionPricer::computeDomainCN() {
	initDomain();
	double w = 1.2;
	double tol = pow(10, -6);

	for (int i = 0; i < N-1; ++i) {
		A[i][i] = 1.0 + alpha;
	}
	for (int i = 1; i < N-2; ++i) {
		A[i][i+1] = -0.5 * alpha;
		A[i][i-1] = -0.5 * alpha;
		A[i-1][i] = -0.5 * alpha;
	}

	// update domain_matrix
	for (int m = 1; m < M + 1; ++m) {
		vector<double> b_vec(N-1);
        for (int n = 1; n < N - 2; ++n) {
            b_vec[n] = 0.5 * alpha * domain[m - 1][n + 2]
            	+ (1 - alpha) * domain[m - 1][n + 1]
            	+ 0.5 * alpha * domain[m - 1][n];
        }
        b_vec[0] = 0.5*alpha * domain[m - 1][2]
        		+ (1 - alpha) * domain[m - 1][1]
        		+ 0.5 * alpha * domain[m - 1][0]
        		+ 0.5 * alpha * domain[m][0];
        b_vec[N - 2] = 0.5 * alpha * domain[m - 1][N]
        		+ (1 - alpha)*domain[m - 1][N - 1]
        		+ 0.5 * alpha * domain[m - 1][N - 2]
        		+ 0.5 * alpha * domain[m][N];
        
		vector<double> x(N + 1);
		x[0] = 0;
		x[N] = 0;
		for (int i = 1; i <= N - 1; ++i) {
			x[i] = option_base.S - option_base.K;
		}

		// find domain[m][:]
		double error = std::numeric_limits<double>::max();
		while (error > tol) {
			vector<double> x_tmp(N + 1);
			x_tmp[0] = 0;
			x_tmp[N] = 0;
			for (int j = 1; j <= N - 1; ++j) {
				x_tmp[j] = (1 - w) * x[j]
							+ 0.5 * w * alpha * (x_tmp[j - 1] + x[j + 1])
							/ (1 + alpha) + w * b_vec[j - 1] / (1 + alpha);
			}
			error = 0.0;
			for (int i = 1; i <= N - 1; ++i) {
				error += pow(x_tmp[i] - x[i], 2);
			}
			error = sqrt(error);
			x = x_tmp;
		}
		// update domain[m][:]
		for (int j = 1; j <= N-1; ++j) { domain[m][j] = x[j]; }
	}
	return;
}

// ---------------------------------------------------
// Compute option stats
// ---------------------------------------------------
void FiniteDifferenceDaoOptionPricer::computeOption(Method method) {
	a = (option_base.r - option_base.q) / option_base.v / option_base.v - 0.5;
	b = pow(a + 1, 2.0) + 2.0 * option_base.q / (option_base.v * option_base.v);

	if (method == FE) {
		computeDomainFE();
		computeOption();
	} else if (method == BE) {
		computeDomainBE();
		computeOption();
	} else if (method == CN) {
		computeDomainCN();
		computeOption();
	}
}

void FiniteDifferenceDaoOptionPricer::computeOption() {
	double exact = BlackScholesPrice(option_base);

	double u = domain[M][N_left];
	double V_approx = exp(-a * x_compute - b * tau_final) * u;
	double error_pointwise = abs(exact - V_approx);

	double x_up = x_compute + delta_x;
	double x_mid = x_compute;
	double x_down = x_compute - delta_x;

	double S_up = option_base.K * exp(x_up);
	double S_mid = option_base.K * exp(x_mid);
	double S_down = option_base.K * exp(x_down);

	double V_up = exp(-a*x_up - b*tau_final) * domain[M][N_left + 1];
	double V_mid = exp(-a*x_mid - b*tau_final) * domain[M][N_left];
	double V_down = exp(-a*x_down - b*tau_final) * domain[M][N_left - 1];
	
	double delta = (V_up - V_down) / (S_up - S_down);
	double gamma = 2 * ((S_mid - S_down) * V_up - (S_up - S_down) * V_mid
					+ (S_up - S_mid) * V_down) / ((S_mid - S_down) * (S_up - S_mid) * (S_up - S_down));
	double delta_t = 2 * delta_tau / option_base.v / option_base.v;
	double V_approx_2 = exp(-a * x_mid - b * (tau_final - delta_tau)) * domain[M-1][N_left];
	// double theta = (V_approx - V_approx_2) / delta_t;
	double theta = (V_approx_2 - V_approx) / delta_t;

	option_stats = vector<double>{ u , V_approx, error_pointwise, delta, gamma, theta};
}

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------
double FiniteDifferenceDaoOptionPricer::f(double x) {
	return option_base.K * exp(a * x) * max(exp(x) - 1, (double)0.0);
}

double FiniteDifferenceDaoOptionPricer::g_right(double tau) {
	return option_base.K * exp(a * x_right + b * tau)
			* (exp(x_right - 2 * option_base.q * tau / pow(option_base.v, 2))
		    - exp(-2 * option_base.r * tau / pow(option_base.v, 2)));
}

