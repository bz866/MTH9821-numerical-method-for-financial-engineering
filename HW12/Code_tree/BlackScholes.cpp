#include "BlackScholes.hpp"
#include "OptionBase.hpp" 
#include <cmath>
#include <math.h>
#include <vector>

// PDF
double n(double x) {
    return 1.0 / sqrt(2.0 * M_PI) * exp(-0.5 * x * x);
}

// CDF
double N(double t) {
    return 0.5 * erfc(-t / sqrt(2));
}

// Option's Black-Scholes price
double BlackScholesPrice(Option opt) {
    return BlackScholesPrice(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type);
}


// Option's Black-Scholes price
double BlackScholesPrice(BarrierOption opt) {
    double a = (opt.r - opt.q) / pow(opt.v, 2) - 0.5;
    double C1 = BlackScholesPrice(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type);
    double C2 = BlackScholesPrice(pow(opt.B,2) / opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type);
    return C1 - pow(opt.B / opt.S, 2 * a) * C2;
}

double BlackScholesPrice(double S, double K, double T, double v, double q, double r, char type) {
    double d1 = (log(S/K) + (r - q + 0.5 * (v * v)) * T)/(v * sqrt(T));
    double d2 = d1 - v * sqrt(T);
    return (type=='C') ? ((S * exp(-q * T) * N(d1)) - (K * exp(-r * T) * N(d2))) :
        ((K * exp(-r * T) * N(-d2))-(S * exp(-q * T) * N(-d1)));
}


double BlackScholesPriceDao(double S, double K, double B, double T, double v, double q, double r, char type) {
    double a = (r - q) / pow(v, 2) - 0.5;
    double C1 = BlackScholesPrice(S, K, T, v, q, r, type);
    double C2 = BlackScholesPrice(pow(B,2) / S, K, T, v, q, r, type);
    return C1 - pow(B / S, 2 * a) * C2;
}

double Delta(double T, double K, double v, double r, double q, double S, char type) {
    double d1 = (log(S/K) + (r - q + 0.5 * (v * v)) * T)/(v * sqrt(T));
    return (type == 'C') ? exp(-q*T) * N(d1) : -exp(-q*T) * N(-d1);
}

double Gamma(double T, double K, double v, double r, double q, double S) {
    double d1 = (log(S/K) + (r - q + 0.5 * (v * v)) * T)/(v * sqrt(T));
    return (n(d1) * exp(-q*T)) / (S * v * sqrt(T));
}

double Theta(double T, double K, double v, double r, double q, double S, char type) {
    double d1 = (log(S/K) + (r - q + 0.5 * (v * v)) * T)/(v * sqrt(T));
    double d2 = d1 - v * sqrt(T);
    return (type == 'C') ? (-n(d1) * S * v * exp(-q * T)/(2 * sqrt(T)) + q * S * exp(-q * T) * N(d1) - r * K * exp(-r * T) * N(d2)) : 
            -n(d1) * S * v * exp(-q * T)/(2 * sqrt(T)) - q * S * exp(-q * T) * N(-d1) + r * K * exp(-r * T) * N(-d2);
}
double Vega(double T, double K, double v, double r, double q, double S)
{
    double d1 = (log(S/K) + (r - q + 0.5 * (v * v)) * T)/(v * sqrt(T));
    return (S * n(d1) * exp(-q*T) * sqrt(T));
}

double BlackScholesGreeks(double S, double K, double T, double v, double q, double r, char type, char Greek) {
    double greek = std::numeric_limits<double>::quiet_NaN();
    switch(Greek) {
        case('D'): greek = Delta(T, K, v, r, q, S, type); break;
        case('G'): greek = Gamma(T, K, v, r, q, S); break;
        case('T'): greek = Theta(T, K, v, r, q, S, type); break;
        case('V'): greek = Vega(T, K, v, r, q, S); break;
    }
    return greek;
}

double BlackScholesImpliedVol(BarrierOption opt, double V, double threshold) {
    double impVol0 = opt.v;
    double impVol = impVol0 + 2 * threshold;
    double BSprice, vega;
    while (abs(impVol - impVol0) > threshold) {
        impVol = impVol0;
        BSprice = BlackScholesPrice(opt.S, opt.K, opt.T, impVol0, opt.q, opt.r, opt.type);
        vega = Vega(opt.T, opt.K, impVol0, opt.r, opt.q, opt.S);
        impVol0 = impVol0 - (BSprice - V) / vega;
    }
    return impVol0;
}

std::vector<double> PriceAndGreeks(BarrierOption opt) {
    double V = BlackScholesPrice(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type);
    double delta = BlackScholesGreeks(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type, 'D');
    double gamma = BlackScholesGreeks(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type, 'G');
    double theta = BlackScholesGreeks(opt.S, opt.K, opt.T, opt.v, opt.q, opt.r, opt.type, 'T');
    return std::vector<double>{V, delta, gamma, theta};
}