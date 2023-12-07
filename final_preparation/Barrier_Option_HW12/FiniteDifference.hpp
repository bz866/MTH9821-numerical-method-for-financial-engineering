#ifndef FINITE_DIFFERENCE_HPP
#define FINITE_DIFFERENCE_HPP

#include <vector>
#include <string>
#include <iostream>
#include "OptionBase.hpp"
#include "BlackScholes.hpp"
#include "LinAlg.h"

enum Method {
	FE = 1,  // Forward Euler
	BE = 2,  // Backward Euler
	CN = 3,  // Crank Nicolson with SOR
};

class FiniteDifferenceDaoOptionPricer {
private:
	void computeDomainFE();
	void computeDomainBE();
	void computeDomainCN();
	void computeOption();
	void initDomain();

	// Domain discretization
	void discretize(double alpha_tmp);

	// Boundary conditions
	double f(double x);
	double g_right(double tau);

	double N_left;
	double N_right;
	double a;
	double b;
	std::vector<double> x_pos;
	std::vector<double> tau_pos;
	BarrierOption option_base;

public:
	FiniteDifferenceDaoOptionPricer(
		BarrierOption option_base_,
		double alpha_tmp_,
		int M_) : option_base(option_base_), M(M_) {
		discretize(alpha_tmp_);
	}

	void computeOption(Method method);

	int M;
	int N;
	double alpha;
	double x_compute;
	double x_left;
	double x_right;
	double delta_tau;
	double delta_x;
	double tau_final;
	std::vector<double> option_stats;

	// Forward Euler
	std::vector<std::vector<double>> domain;

	// Backward Euler
	vector<vector<double>> A;

	// Crank-Nicolson
	vector<vector<double>> B;
};


#endif /* FINITE_DIFFERENCE_HPP */
