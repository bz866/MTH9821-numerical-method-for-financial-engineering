#ifndef BS_HPP
#define BS_HPP

#include <string>
#include <vector>

// black scholes formula to calcuate call and put price and greeks 
std::vector<double> BlackScholes(double S, double K,
	double r, double T, double q, double sig, const char option);

// black scholes formula to calcuate down and out call value 
double BlackScholes_DaO(double S, double K,
	double r, double T, double q, double sig, double B, const char option);

#endif // !BS_HPP
