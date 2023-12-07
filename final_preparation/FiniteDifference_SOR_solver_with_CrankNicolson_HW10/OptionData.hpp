#ifndef OptionData_hpp
#define OptionData_hpp

#include <iostream>

struct option
{
    double S0;     // Spot price
    double K;      // Strike price
    double T;      // Maturity
    double v;      // Volatility
    double q;      // Dividends
    double r;      // Risk-free interest rate
    char AmerEuro; // 'A' for American and 'E' for European options
    char PutCall;  // 'P' for put and 'C' for call options
};

#endif /* OptionData_hpp */
