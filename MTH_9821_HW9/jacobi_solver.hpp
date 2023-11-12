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

#endif