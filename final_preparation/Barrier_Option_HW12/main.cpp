#include "FiniteDifference.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "LinAlg.h"
#include "BS.hpp"
#include "Tree.hpp"
#include "BlackScholes.hpp"
#include "PathDependentOption.hpp"

using namespace std;

int main() 
{
    cout << "Question 1: Trinomial Trees pricing" << endl;
    {
        BarrierOption option_base = {
            .S = 42.0,
            .K = 40.0,
            .B = 36.0,
            .T = 7.0 / 12,
            .v = 0.28,
            .q = 0.015,
            .r = 0.04,
            .type = 'C'
        };

        ofstream txt("Q1_Trinomial_Trees.txt", std::ofstream::out);
        if (txt.is_open())
        {
            cout << setprecision(10);
            cout << "question 1" << "\n";
            txt << setprecision(10);
            txt << "question 1" << "\n";
            auto bsvalue = BlackScholes_DaO(option_base.S, option_base.K, option_base.r, option_base.T, option_base.q, option_base.v, option_base.B, option_base.type);
            cout << "part 1" << endl;
            txt << "part 1" << "\n";
            cout << "bsvalue:" << bsvalue << endl;
            txt << "bsvalue:" << bsvalue << "\n";
            cout << endl;
            txt << "\n";

            cout << "part 2" << endl;
            txt << "part 2" << "\n";
            cout << "print out approxmiate errors in Q1_Trinomial_Trees.txt, plot in excel" << endl;
            for (int N = 10; N < 1001; N++)
            {
                auto result = TT_euro_DAOC(N, option_base.S, option_base.K, option_base.r, option_base.T, option_base.q, option_base.v, option_base.B, option_base.type);
                txt << N << "\t" << abs(result[0] - bsvalue) << "\n";
            }
            cout << endl;
            txt << "\n";

            cout << "part 3" << endl;
            txt << "part 3" << "\n";

            vector<double> k_list;
            vector<double> N_k_list;
            for (int k = 2; k < 31; k++)
            {
                int N_k_temp = floor(k * k * 3 * option_base.v * option_base.v * option_base.T / log(option_base.S / option_base.B) / log(option_base.S / option_base.B));
                auto result = TT_euro_DAOC(N_k_temp, option_base.S, option_base.K, option_base.r, option_base.T, option_base.q, option_base.v, option_base.B, option_base.type);

                cout << k << "\t" << N_k_temp << "\t" << abs(result[0] - bsvalue) << endl;
                txt << k << "\t" << N_k_temp << "\t" << abs(result[0] - bsvalue) << endl;


            }
        }
    }

    cout << "Question 2: Monte Carlo pricing (in Q2_monte_carlo.txt file)" << endl;
    {
        BarrierOption option_base = {
            .S = 42.0,
            .K = 40.0,
            .B = 36.0,
            .T = 7.0 / 12,
            .v = 0.28,
            .q = 0.015,
            .r = 0.04,
            .type = 'C'
        };

        unsigned long int m = 200, n;
        ofstream txt("Q2_monte_carlo.txt", std::ofstream::out);
        if (txt.is_open())
        {
            txt << "Monte Carlo method" << "\n";
            txt << setprecision(10);
            double N_k;
            // double V;
            PathDependentOption barrier_option = PathDependentOption(option_base);
            double CF_Cdao = barrier_option.ClosedFormCdao();
            double MC_Cdao;

            txt << "-----------[ Part 1 ]-----------\n"
                << "m=200, n, V(n), |CF_Cdao-V(n)|"
                << endl;
            for (int k = 0; k < 10; ++k) {
                N_k = 10000 * pow(2,k);
                MC_Cdao = barrier_option.MonteCarloCdao(N_k / m, m);
                txt << m << "\t" << N_k / m << "\t" << MC_Cdao << "\t"
                         << std::abs(CF_Cdao - MC_Cdao) << endl;
            }
            txt << "\n-----------[ Part 2 ]-----------\n"
                << "m_k, n_k, V(n), |CF_Cdao-V(n_n)|"
                << endl;
            for (int k = 0; k < 10; ++k) {
                N_k = 10000 * pow(2, k);
                m = (unsigned long int)ceil(pow((double)N_k, 1./3.) * pow(option_base.T, 2./3.));
                n = (unsigned long int)floor(N_k / m);
                MC_Cdao = barrier_option.MonteCarloCdao(n, m);
                
                txt << m <<  "\t" << n << "\t" << MC_Cdao << "\t"
                         << std::abs(CF_Cdao - MC_Cdao) << endl;
            }
        }
    }

    cout << "Question 3: Finite differences (in Q3_FD.txt file)" << endl;
    {
        BarrierOption option_base = {
            .S = 42.0,
            .K = 40.0,
            .B = 36.0,
            .T = 7.0 / 12,
            .v = 0.28,
            .q = 0.015,
            .r = 0.04,
            .type = 'C'
        };

        ofstream txt("Q3_FD.txt", std::ofstream::out);
        if (txt.is_open())
        {
            double alpha_tmp = 0.4;
            int M;
            txt << setprecision(10);
            txt << "Domain discretization with α_tmp = 0.4" << "\n";
            txt << "M,   alpha,   x_left,  x_right, N,   δx,  δτ\n";
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                txt << M << ",   " << pricer.alpha << ",   " << pricer.x_left << ",  "
                    << pricer.x_right << "\t" << pricer.N << ",   " << pricer.delta_x << "\t"
                    << pricer.delta_tau << endl;
            }

            txt << endl;
            txt << "Domain discretization with α_tmp = 4" << "\n";
            txt << "M,   alpha,   x_left,  x_right, N,   δx,  δτ\n";
            alpha_tmp = 4.0;
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                txt << M << ",   " << pricer.alpha << ",   " << pricer.x_left << ",  "
                    << pricer.x_right << "\t" << pricer.N << ",   " << pricer.delta_x << "\t"
                    << pricer.delta_tau << endl;
            }

            txt << endl;
            txt << "Forward Euler with α_tmp = 0.4" << "\n";
            txt << "M,   u value, Option value,    Pointwise Error, Delta_central,   Gamma_central,   Theta_forward\n";
            alpha_tmp = 0.4;
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::FE);
                txt << M << ",   " << pricer.option_stats[0] << "\t" << pricer.option_stats[1] << "\t" << pricer.option_stats[2] << "\t"
                << pricer.option_stats[3] << "\t" << pricer.option_stats[4] << "\t" << pricer.option_stats[5] << endl;
            }

            txt << endl;
            txt << "Backward Euler with LU and α_tmp = 0.4" << "\n";
            txt << "M,   u value, Option value,    Pointwise Error, Delta_central,   Gamma_central,   Theta_forward\n";
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::BE);
                txt << M << ",   " << pricer.option_stats[0] << "\t" << pricer.option_stats[1] << "\t" << pricer.option_stats[2] << "\t"
                << pricer.option_stats[3] << "\t" << pricer.option_stats[4] << "\t" << pricer.option_stats[5] << endl;
            }

            txt << endl;
            txt << "Backward Euler with LU and α_tmp = 4" << "\n";
            txt << "M,   u value, Option value,    Pointwise Error, Delta_central,   Gamma_central,   Theta_forward\n";
            alpha_tmp = 4;
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::BE);
                txt << M << ",   " << pricer.option_stats[0] << "\t" << pricer.option_stats[1] << "\t" << pricer.option_stats[2] << "\t"
                << pricer.option_stats[3] << "\t" << pricer.option_stats[4] << "\t" << pricer.option_stats[5] << endl;
            }

            txt << endl;
            txt << "Crank-Nicolson with SOR and α_tmp = 0.4" << "\n";
            txt << "M,   u value, Option value,    Pointwise Error, Delta_central,   Gamma_central,   Theta_forward\n";
            alpha_tmp = 0.4;
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::CN);
                txt << M << ",   " << pricer.option_stats[0] << "\t" << pricer.option_stats[1] << "\t" << pricer.option_stats[2] << "\t"
                << pricer.option_stats[3] << "\t" << pricer.option_stats[4] << "\t" << pricer.option_stats[5] << endl;
            }

            txt << endl;
            txt << "Crank-Nicolson with SOR and α_tmp = 4" << "\n";
            txt << "M,   u value, Option value,    Pointwise Error, Delta_central,   Gamma_central,   Theta_forward\n";
            alpha_tmp = 4;
            for (int i = 1; i < 5; ++i) {
                M = pow(4, i);
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::CN);
                txt << M << ",   " << pricer.option_stats[0] << "\t" << pricer.option_stats[1] << "\t" << pricer.option_stats[2] << "\t"
                << pricer.option_stats[3] << "\t" << pricer.option_stats[4] << "\t" << pricer.option_stats[5] << endl;
            } 

            txt << endl;
            txt << "Forward Euler with α_tmp = 0.4" << "\n";
            txt << "m,   x_left,  x_compute,   x_1, x_2, x_3, x_right\n";
            alpha_tmp = 0.4;
            M = 4;
            for (int m = 0; m <= M; ++m) {
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::FE);
                txt << m << "\t" ;
                for (auto iter = pricer.domain[m].begin(); iter != pricer.domain[m].end(); ++iter) {
                    txt << *iter << "\t" ;
                }
                txt << endl;
            }

            txt << endl;
            txt << "Backward Euler with LU and α_tmp = 0.4" << "\n";
            txt << "m,   x_left,  x_compute,   x_1, x_2, x_3, x_right\n";
            M = 4;
            for (int m = 0; m <= M; ++m) {
                FiniteDifferenceDaoOptionPricer pricer(option_base, alpha_tmp, M);
                pricer.computeOption(Method::BE);
                txt << m << "\t" ;
                for (auto iter = pricer.domain[m].begin(); iter != pricer.domain[m].end(); ++iter) {
                    txt << *iter << "\t" ;
                }
                txt << endl;
            }
        }
    }
    

    return 0;
}
