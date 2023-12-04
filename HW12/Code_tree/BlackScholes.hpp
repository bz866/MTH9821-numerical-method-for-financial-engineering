#ifndef BLACK_SCHOLES_HPP
#define BLACK_SCHOLES_HPP

#include "OptionBase.hpp"
#include <vector>
#include <iostream> 


double BlackScholesPrice(Option opt);

double BlackScholesPrice(BarrierOption opt);

double BlackScholesPrice(double S, double K, double T, double v, double q, double r, char type);

double BlackScholesPriceDao(double S, double K, double B, double T, double v, double q, double r, char type);

double BlackScholesGreeks(double S, double K, double T, double v, double q, double r, char type, char Greek);

double BlackScholesImpliedVol(BarrierOption opt, double V, double threshold);

double Delta(double T, double K, double v, double r, double q, double S, char type);

double Gamma(double T, double K, double v, double r, double q, double S);

double Theta(double T, double K, double v, double r, double q, double S, char type);

double Vega(double T, double K, double v, double r, double q, double S);

std::vector<double> PriceAndGreeks(BarrierOption opt);

#endif /* BLACK_SCHOLES_HPP */
