#ifndef OPTION_BASE_HPP
#define OPTION_BASE_HPP

#include <iostream>


// Down–and–Out barrier option

struct Option 
{
    double S;  // Spot price
    double K;  // Strike price
    double T;  // Expiry
    double v;  // Volatility
    double q;  // Dividend yield
    double r;  // Risk-free rate
    char type; // 'C' for Call, 'P' for Put
};

struct BarrierOption
{
    double S;  // Spot price
    double K;  // Strike price
    double B;  // Barrier
    double T;  // Expiry
    double v;  // Volatility
    double q;  // Dividend yield
    double r;  // Risk-free rate
    char type; // 'C' for Call, 'P' for Put
};

#endif /* OPTION_BASE_HPP */
