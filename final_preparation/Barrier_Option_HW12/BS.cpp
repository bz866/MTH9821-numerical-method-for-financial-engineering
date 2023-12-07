#include "BS.hpp"
#include <cmath>

// black scholes formula to calcuate call and put price and greeks 
std::vector<double> BlackScholes(double S, double K,
	double r, double T, double q, double sig, const char option)
{
	double d1 = (std::log(S / K) + (r - q + sig * sig / 2) * T) / (sig * std::sqrt(T));
	double d2 = d1 - sig * std::sqrt(T);

	double pi = 3.141592653589793238462;
	double Nd1 = 0.5 * std::erfc(-d1 * std::sqrt(.5f));
	double Nd2 = 0.5 * std::erfc(-d2 * std::sqrt(.5f));

	double price = 0;
	double delta = 0;
	double vega = S * std::exp(-q * T) / std::sqrt(2.0 * pi) * std::exp(-0.5 * d1 * d1) * std::sqrt(T);
	double gamma = std::exp(-q * T) / (S * sig * std::sqrt(T)) / std::sqrt(2.0 * pi) * std::exp(-0.5 * d1 * d1);
	double theta = 0;
	if (option == 'C')
	{
		price = Nd1 * S * std::exp(-q * T) - Nd2 * K * std::exp(-r * (T));
		delta = std::exp(-q * T) * Nd1;
		//theta= -(S * sig * exp(-q * T)) / (2 * sqrt(2 * pi * T)) * exp(-d1 * d1 / 2.0) - 
		//	q * S * exp(-q * T) * Nd1 + r * K * exp(-r * T) * Nd2;
		theta= -(S * sig * exp(-q * T)) / (2 * sqrt(2 * pi * T)) * exp(-d1 * d1 / 2.0) + 
			q * S * exp(-q * T) * Nd1 - r * K * exp(-r * T) * Nd2;
	}
	if (option == 'P')
	{
		price = (1 - Nd2) * K * std::exp(-r * (T)) - (1 - Nd1) * S * std::exp(-q * T);
		delta = std::exp(-q * T) * (Nd1 - 1);
		theta = -(S * sig * exp(-q * T)) / (2 * sqrt(2 * pi * T)) * exp(-d1 * d1 / 2.0) -
			q * S * exp(-q * T) * (1-Nd1) + r * K * exp(-r * T) * (1-Nd2);
	}
	
	return std::vector<double>{price, delta, vega, gamma, theta};

	
	
}

// black scholes formula to calcuate down and out call value 
double BlackScholes_DaO(double S, double K,
	double r, double T, double q, double sig, double B, const char option)
{
	double a = (r - q) / sig / sig - 0.5;
	double call_value = BlackScholes(S, K, r, T, q, sig, option)[0];
	double barrier_call = BlackScholes(B * B / S, K, r, T, q, sig, option)[0];
	double result = call_value - std::pow((B / S), 2 * a) * barrier_call;
	return result;
}

