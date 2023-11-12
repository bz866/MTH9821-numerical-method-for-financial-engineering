#ifndef SOR_SOLVER_HPP
#define SOR_SOLVER_HPP

#include "stopping_criteria.hpp"
#include <Dense>
#include <tuple>
#include <cmath>

// SOR iterative solver function
std::tuple<vec, unsigned int> sor(const mat& A, const vec& b, double omega, double tolerance, StoppingCriterion criterion) {
    vec x = vec::Zero(A.rows());  // Initial guess x_0 as zero vector
    vec x_prev = x;
    unsigned int iterations = 0;

    do {
        x_prev = x;
        for (int i = 0; i < A.rows(); ++i) {
            double sigma1 = 0.0;
            for (int j = 0; j < i; ++j) {
                sigma1 += A(i, j) * x[j];
            }

            double sigma2 = 0.0;
            for (int j = i + 1; j < A.cols(); ++j) {
                sigma2 += A(i, j) * x_prev[j];
            }

            x[i] = (1 - omega) * x_prev[i] + (omega / A(i, i)) * (b[i] - sigma1 - sigma2);
        }

        if (criterion == consecutive && isConsecutiveClose(x, x_prev, tolerance)) {
            break;
        }
        else if (criterion == residual && isResidualSmall(A, b, x, tolerance)) {
            break;
        }

        iterations++;
    } while (true);

    return std::make_tuple(x, iterations);
}

#endif