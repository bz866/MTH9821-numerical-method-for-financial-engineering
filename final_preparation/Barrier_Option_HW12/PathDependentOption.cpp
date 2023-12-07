#include <cmath>
#include <vector>
#include "PathDependentOption.hpp"
#include "BlackScholes.hpp"
#include "NumericalMethod.hpp"

double PathDependentOption::ClosedFormCdao() {
    return BlackScholesPrice(option_base);
    // double a = (option_base.r - option_base.q) / (option_base.v * option_base.v) - 0.5;
    // double barrier_call = BlackScholesPrice(option_base.B * option_base.B / option_base.S, option_base.K, option_base.T, option_base.v, option_base.q, option_base.r, 'C');
    // return BlackScholesPrice(option_base) - pow(option_base.B / option_base.S, 2 * a) * barrier_call;
}

double PathDependentOption::MonteCarloCdao(unsigned long int n, unsigned long int m) {
    /* 
        n = the number of paths
        m = the number of time steps on each path
    */
    double result = 0;
    double dt = option_base.T / double(m);

    NumericalMethod numericalMethod;
    std::vector<long double> randomVec = numericalMethod.InverseTransformRandom(n*m);

    for (unsigned long int i = 0; i < n; ++i) {
        double S = option_base.S;
        for (unsigned long int j = 0; j < m; ++j) {
            S = S * exp((option_base.r - option_base.q - 0.5 * option_base.v * option_base.v) * dt +
                option_base.v * sqrt(dt) * randomVec[i * m + j]);
            if (S <= option_base.B) {
                S = 0;
                break;
            }
        }
        if (S > option_base.K) {
            result += exp(-option_base.r * option_base.T) * (S - option_base.K);
        }
    }

    return result / n;
}
