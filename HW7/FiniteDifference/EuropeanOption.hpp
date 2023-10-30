#ifndef EUROPEAN_OPTION_HPP
#define EUROPEAN_OPTION_HPP

#include <cmath>
#include <numbers>
#include <iostream>
#include <algorithm>

class EuropeanOption 
{
private:
    double t_; // current time
    double S_; // spot price
    double K_; // strike price
    double T_; // maturity
    double sigma_; // volatility
    double r_; // const interest rate
    double q_; // dividend rate

    // Black-Scholes pricing values
    double d1_;
    double d2_;
    double Zd1_;
    double Zd2_;
    double Nd1_;
    double Nd2_;
    double q_disc_; // Discount by dividend rate: exp(-q(T-t))
    double r_disc_; // Discount by interest rate: exp(-r(T-t))

    double Phi(double t) const
    {
        return std::erfc(-t / std::sqrt(2.)) / 2.;
    }

    double Z(double t) const
    {
        return std::exp(-t * t / 2.) / std::sqrt(2. * M_PI);
    }

public:
    EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q);
    ~EuropeanOption() = default;


    // properties defined inline
    double Call() const
    {
        return S_ * q_disc_ * Nd1_ - K_ * r_disc_ * Nd2_;
    }

    double Put() const
    {
        return -S_ * q_disc_ * (1. - Nd1_) + K_ * r_disc_ * (1. - Nd2_);
    }

    double CallPayoff(const double S, const double t) const
    {
        return std::max(S - K_, 0.);
    }

    double PutPayoff(const double S, const double t) const
    {
        return std::max(K_ - S, 0.);
    }

    // Greeks
    double DeltaCall() const
    {
        return q_disc_ * Nd1_;
    }
    double DeltaPut() const
    {
        return -q_disc_ * (1. - Nd1_);
    }

    double GammaCall() const
    {
        return q_disc_ / (S_ * sigma_ * std::sqrt(T_ - t_)) * Zd1_;
    }
    double GammaPut() const
    {
        return this->GammaCall();
    }

    double ThetaCall() const
    {
        double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;

        res += q_ * S_ * q_disc_ * Nd1_;

        res += -r_ * K_ * r_disc_ * Nd2_;

        return res;
    }

    double ThetaPut() const
    {
        double res = -(S_ * sigma_ * q_disc_) / (2. * std::sqrt(T_ - t_)) * Zd1_;

        res += -q_ * S_ * q_disc_ * (1 - Nd1_);

        res += r_ * K_ * r_disc_ * (1 - Nd2_);

        return res;
    }

    double VegaCall() const
    {
        return S_ * q_disc_ * std::sqrt(T_ - t_) * Zd1_;
    }

    double VegaPut() const
    {
        return this->VegaCall();
    }

    double RhoCall() const
    {
        return K_ * (T_ - t_) * r_disc_ * Nd2_;
    }

    double RhoPut() const
    {
        return -K_ * (T_ - t_) * r_disc_ * (1. - Nd2_);
    }
};



#endif
