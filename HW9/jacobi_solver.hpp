#ifndef JACOBI_SOLVER_HPP
#define JACOBI_SOLVER_HPP

#include "stopping_criteria.hpp"
#include <Dense>
#include <tuple>
#include <cmath>

// Assuming the provided enum, type definitions, and helper functions

// Jacobi iterative solver function
std::tuple<vec, unsigned int> jacobi(const mat& A, const vec& b, const vec& x_0, double tolerance, StoppingCriterion criterion) {
    vec x = x_0;
    vec new_x = x;
    unsigned int iterations = 0;

    do {
        for (int i = 0; i < A.rows(); ++i) {
            double sigma = 0.0;
            for (int j = 0; j < A.cols(); ++j) {
                if (i != j) {
                    sigma += A(i, j) * x[j];
                }
            }
            new_x[i] = (b[i] - sigma) / A(i, i);
        }

        if (criterion == consecutive && isConsecutiveClose(new_x, x, tolerance)) {
            break;
        }
        else if (criterion == residual && isResidualSmall(A, b, new_x, tolerance)) {
            break;
        }

        x = new_x;
        iterations++;
    } while (true);

    return std::make_tuple(x, iterations);
}

// Assuming the provided enum, type definitions, and helper functions

// Modified jacobiBanded function to include stopping criteria
std::tuple<vec, unsigned int> jacobiBanded(const mat& A, const vec& b, int bandwidth, double tolerance, StoppingCriterion criterion) {
    int n = A.rows();
    vec x = vec::Zero(n);
    vec new_x = x;
    vec x_prev = x;
    unsigned int iterations = 0;

    while (true) {
        x_prev = x;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = std::max(0, i - bandwidth); j <= std::min(n - 1, i + bandwidth); ++j) {
                if (j != i) {
                    sum += A(i, j) * x[j];
                }
            }
            new_x[i] = (b[i] - sum) / A(i, i);
        }

        // Check for convergence
        if ((criterion == consecutive && isConsecutiveClose(new_x, x, tolerance)) ||
            (criterion == residual && isResidualSmall(A, b, new_x, tolerance))) {
            break;
        }

        x = new_x;
        iterations++;
    }

    return std::make_tuple(x, iterations);
}

#endif