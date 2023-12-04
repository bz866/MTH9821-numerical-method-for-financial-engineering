#include <cmath>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include "BS.hpp"
#include "Tree.hpp"
#include "BlackScholes.hpp" 

using namespace std;

int main() 
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
    cout << "Question 1: Trinomial Trees pricing" << endl;
    {
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

    return 0;
}
