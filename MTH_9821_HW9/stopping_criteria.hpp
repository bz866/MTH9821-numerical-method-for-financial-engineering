#ifndef STOPPING_CRITERIA_HPP
#define STOPPING_CRITERIA_HPP

#include <vector>
#include <cmath>

enum StoppingCriterion { consecutive, residual };

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1 , -1 , unsigned long> permutation;

// Function to check if consecutive iterations are close in norm
bool isConsecutiveClose(const vec& x, const vec& x_prev, double tolerance) {
    double norm = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        norm += (x[i] - x_prev[i]) * (x[i] - x_prev[i]);
    }
    return std::sqrt(norm) < tolerance;
}

// Function to check if the residual is small
bool isResidualSmall(const mat& A, const vec& b, const vec& x, double tolerance) {
    vec residual = b - A * x; // Compute the residual
    return residual.norm() < tolerance; // Check if the norm is less than the tolerance
}

#endif