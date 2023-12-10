# include <iostream>
# include <math.h>
# include <vector>

using namespace std;


vector <double> EuropeanBinomialPutPercentageDividend(int N, double S0, double K, double r, double g, double dividendtime, double sigma, double T) {

    vector<vector <double>> S(N + 1);
    vector<vector <double>> V(N + 1);

    int size = S.size();

    for (int i = 0; i < size; i++) {

        S[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    for (int j = N; j >= 0;j--){

        S[N][j]= S0 * pow(1 - g, static_cast<int>(N * deltaT / dividendtime))* pow(u, 2 * j - N); 
        V[N][j] = max(K- S[N][j], 0.0);

    }

    for (int i = N-1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S0 * pow(1 - g, static_cast<int>(i * deltaT / dividendtime))* pow(u, 2 * j - i);
            V[i][j] = DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]);

        }
    }

    double value = V[0][0];
    double delta = (V[1][1] - V[1][0]) / (S[1][1] - S[1][0]);
    double gamma = ((V[2][2] - V[2][1]) / (S[2][2] - S[2][1]) - (V[2][1] - V[2][0]) / (S[2][1] - S[2][0])) / ((S[2][2] - S[2][0]) / 2);
    double theta = (V[2][1] - V[0][0]) / (2 * deltaT);

    return { value, delta, gamma, theta };

}

vector <double> AmericanBinomialPutPercentageDividend(int N, double S0, double K, double r, double g, double dividendtime, double sigma, double T) {

    vector<vector <double>> S(N + 1);
    vector<vector <double>> V(N + 1);

    int size = S.size();

    for (int i = 0; i < size; i++) {

        S[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    for (int j = N; j >= 0; j--) {

        S[N][j] = S0 * pow(1 - g, static_cast<int>(N * deltaT / dividendtime)) * pow(u, 2 * j - N);
        V[N][j] = max(K - S[N][j], 0.0);

    }

    for (int i = N - 1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S0 * pow(1 - g, static_cast<int>(i * deltaT / dividendtime)) * pow(u, 2 * j - i);
            V[i][j] = max(DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]), K - S[i][j]);

        }
    }

    double value = V[0][0];
    double delta = (V[1][1] - V[1][0]) / (S[1][1] - S[1][0]);
    double gamma = ((V[2][2] - V[2][1]) / (S[2][2] - S[2][1]) - (V[2][1] - V[2][0]) / (S[2][1] - S[2][0])) / ((S[2][2] - S[2][0]) / 2);
    double theta = (V[2][1] - V[0][0]) / (2 * deltaT);

    return { value, delta, gamma, theta };

}


int main() {

    double spotPrice = 50.0;
    double strikePrice = 55.55;
    double volatility = 0.25;
    double dividendPercentage = 0.015;
    double dividendPeriod = 2.0 / 12.0;
    double riskFreeRate = 0.04;
    double optionTerm = 3.0 / 12.0;

    int eurIniN = 6;
    int amerIniN = 6;

    double tolerance = pow(10, -4);

    double eurDiff = 1.0;
    double amerDiff = 1.0;
    
    vector <double> opt_eur1;
    vector <double> opt_eur2;

    vector<double> opt_amer1;
    vector<double> opt_amer2;

    
        while (eurDiff >= tolerance) {

            opt_eur1 = EuropeanBinomialPutPercentageDividend(eurIniN, spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendPeriod, volatility, optionTerm);
            opt_eur2 = EuropeanBinomialPutPercentageDividend(eurIniN*2, spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendPeriod, volatility, optionTerm);
            eurDiff = abs(opt_eur1[0] - opt_eur2[0]);
            eurIniN *= 2;
        }
        
        cout << opt_eur1[0] << "\t" << eurIniN << "\t" << opt_eur1[1] << "\t" << opt_eur1[2] << "\t" << opt_eur1[3] << endl;
                
        while (amerDiff >= tolerance) {

            opt_amer1 = AmericanBinomialPutPercentageDividend(amerIniN, spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendPeriod, volatility, optionTerm);
            opt_amer2 = AmericanBinomialPutPercentageDividend(amerIniN * 2, spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendPeriod, volatility, optionTerm);
            amerDiff = abs(opt_amer1[0] - opt_amer2[0]);
            amerIniN *= 2;
        }
        
        cout << opt_amer1[0] << "\t" << amerIniN << "\t" << opt_amer1[1] << "\t" << opt_amer1[2] << "\t" << opt_amer1[3] << endl;
        
}