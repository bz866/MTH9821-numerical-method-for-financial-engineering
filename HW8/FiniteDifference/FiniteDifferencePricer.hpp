#ifndef FINITE_DIFFERENCE_PRICER_HPP
#define FINITE_DIFFERENCE_PRICER_HPP

#include "EuropeanOption.hpp"
#include "utils.hpp"
#include <iostream>
#include <functional>
#include <cmath>
#include <math.h>
#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <iterator> 

using namespace std;

class FiniteDifferencePricer
{
private:
    // heat equation transformation coefficients
    double a_;
    double b_;

    // finite difference hyperparameters 
    double alpha_temp_;
    std::size_t M_;
    bool flag_alpha_temp_defined_ = false;
    bool flag_M_defined_ = false;

    // finite difference parameters
    double tau_final_;
    double x_l_;
    double x_r_;

    // boundary conditions as lambdas
    std::function<double(double)> boundary_tau_0_;
    std::function<double(double)> boundary_x_l_;
    std::function<double(double)> boundary_x_r_;
    
    // private computational utilies
    // European Part
    // V(S, t) = exp(-a*x-b*tau) * u(x, tau)
    std::function<double(double, double, double, double, double)>compute_V_from_heat_equation_solution();
    std::tuple<double, std::size_t, double, double> compute_domain_params(const std::size_t M) const;
    std::tuple<std::vector<double>, double> buildMesh(const std::size_t N, const double d_x);
    std::vector<double> computeNextUMesh(
        const double tau, const double alpha, const std::vector<double> &x_mesh, const std::vector<double> &u_mesh);
    std::pair<std::vector<double>, std::vector<double>> buildUMeshFromFiniteDifference(
        const double alpha, const std::vector<double>& x_mesh, std::size_t M, double d_tau);
    std::vector<double> approximateSecurityValue(
        const std::vector<double> &x_mesh, const std::vector<double>& u_mesh, const double d_x);
    double computeRMSError(
        const std::vector<double>& x_mesh, const std::vector<double>& u_mesh);
    std::vector<double> computeGreeks(
        std::size_t interval_i, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh, 
        std::vector<double>& u_mesh_M_minus_1_th, double d_tau, double V_approx);
    
    // American Part
    std::vector<double> computeNextUMeshAmerican(
        const double tau, const double alpha, const std::vector<double>& x_mesh, const std::vector<double>& u_mesh);
    std::pair<std::vector<double>, std::vector<double>> buildUMeshFromFiniteDifferenceAmericanPut(
        const double alpha, const std::vector<double>& x_mesh, std::size_t M, double d_tau);
    double varianceReductionAmericanPut(double V_approx, std::size_t M);
    // American Early Exercise
    std::vector<double> compute_Sopt(
        const double alpha, const std::vector<double> &x_mesh, std::size_t M, double d_tau);
    std::array<std::vector<double>, 2> AmericanPutEarlyExercise(std::size_t M);

    // European Put with Backward Euler and Crank Nicolson
    std::vector<double> tridiagonalSolve(
        const std::vector<double>& l,
        const std::vector<double>& d,
        const std::vector<double>& u,
        const std::vector<double>& rhs);
    std::vector<double> computeNextUMesh_BE(
        const double tau, const double alpha, const std::vector<double> &x_mesh, const std::vector<double> &u_mesh);
    std::pair<std::vector<double>, std::vector<double>> buildUMeshFromFiniteDifference_BE(
        const double alpha, const std::vector<double>& x_mesh, std::size_t M, double d_tau);



    // boundary conditions setters
    void setUpEuropeanCallBoundaryCondition(const std::size_t M);
    void setUpEuropeanPutBoundaryCondition(const std::size_t M);
    void setUpAmericanCallBoundaryCondition(const std::size_t M);
    void setUpAmericanPutBoundaryCondition(const std::size_t M);

public:
    // Option Properties
    double S0_;
    double K_;
    double T_;
    double sigma_;
    double r_;
    double q_;

    FiniteDifferencePricer(double S0, double K, double T, double sigma, double r, double q);
    ~FiniteDifferencePricer() = default;
    
    // hyperparameters setter
    void setHyperparameters(const double alpha_temp, const double M);
    void validateHyperparameters() const;

    // pricers
    std::vector<double> priceEuropeanPut(const std::size_t M);
    std::vector<double> priceAmericanPut(const std::size_t M);
    std::pair<std::vector<double>, std::vector<double>> priceAmericanPutEarlyExercise(const std::size_t M);
    std::vector<double> priceEuropeanPut_BE(const std::size_t M);
};

#endif
