#ifndef GAUSS_SIEDEL_SOLVER_HPP
#define GAUSS_SIEDEL_SOLVER_HPP

#include "stopping_criteria.hpp"
#include <Dense>
#include <tuple>
#include <cmath>


std::tuple<vec, unsigned int> gauss_seidel(const mat& A, const vec& b, const vec& x_0, double tolerance, StoppingCriterion criterion) {
    vec x = vec::Zero(A.rows()); // Initial guess x_0 as zero vector
    vec x_prev = x;
    unsigned int iterations = 0;

    do {
        x_prev = x;
        for (int i = 0; i < A.rows(); ++i) {
            double sigma = 0.0;
            for (int j = 0; j < i; ++j) {
                sigma += A(i, j) * x[j];
            }
            for (int j = i + 1; j < A.cols(); ++j) {
                sigma += A(i, j) * x_prev[j];
            }
            x[i] = (b[i] - sigma) / A(i, i);
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