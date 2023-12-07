#ifndef BlackScholes_hpp
#define BlackScholes_hpp

#include <iostream>

// Put or Call option value from Black-Scholes
double BlackScholes(double S, double K, double T, double v, double q, double r, char PutCall);

// Greeks
double BlackScholesGreeks(double S, double K, double T, double v, double q, double r, char PutCall, char Greek);

// Delta of call
double CallDelta(double T, double K, double v, double r, double q, double S);

// Delta of put
double PutDelta(double T, double K, double v, double r, double q, double S);

// Gamma
double Gamma(double T, double K, double v, double r, double q, double S);

// Vega
double Vega(double T, double K, double v, double r, double q, double S);

// Theta of call
double CallTheta(double T, double K, double v, double r, double q, double S);

// Theta of put
double PutTheta(double T, double K, double v, double r, double q, double S);

#endif /* BlackScholes_hpp */
