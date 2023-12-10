# include <iostream>
# include <math.h>
# include <vector>

using namespace std;


vector <double> EuropeanBinomialPutFixedDividend(int N, double S0, double K, double r, double div, double dividendtime, double sigma, double T) {

    vector<vector <double>> S(N + 1);
    vector<vector <double>> D(N + 1);
    vector<vector <double>> V(N + 1);

    int size = S.size();

    for (int i = 0; i < size; i++) {

        S[i].resize(i + 1);
        D[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    double S_adj;
    double sumDiv = 0.0;
    
    for (int i = 1; i <= static_cast<int>(T / dividendtime); i++) {

        sumDiv += div * exp(-r * i * dividendtime);

    }
    
    S_adj = S0 - sumDiv;

    
    for (int i = 0; i < N + 1; i++) {

        double s = 0.0;

        for (int k = 1; k <= static_cast<int>(T / dividendtime); k++) {

            if (i * deltaT < k * dividendtime) {

                s += div * exp(-r * (k * dividendtime - i * deltaT));

            }

            for (int j = i; j >= 0; j--) {

                D[i][j] = 0.0 + s;

            }
        }
    }

    for (int j = N; j >= 0; j--) {

        S[N][j] = S_adj * pow(u, 2 * j - N)+D[N][j];
        V[N][j] = max(K - S[N][j], 0.0);

    }

    for (int i = N - 1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S_adj * pow(u, 2 * j - i) + D[i][j];
            V[i][j] = DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]);

        }
    }
    
    double value = V[0][0];
    double delta = (V[1][1] - V[1][0]) / (S[1][1] - S[1][0]);
    double gamma = ((V[2][2] - V[2][1]) / (S[2][2] - S[2][1]) - (V[2][1] - V[2][0]) / (S[2][1] - S[2][0])) / ((S[2][2] - S[2][0]) / 2);
    double theta = (V[2][1] - V[0][0]) / (2 * deltaT);

    return { value, delta, gamma, theta };

}




vector <double> AmericanBinomialPutFixedDividend(int N, double S0, double K, double r, double div, double dividendtime, double sigma, double T) {

    vector<vector <double>> S(N + 1);
    vector<vector <double>> D(N + 1);
    vector<vector <double>> V(N + 1);

    int size = S.size();

    for (int i = 0; i < size; i++) {

        S[i].resize(i + 1);
        D[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    double S_adj;
    double SumDiv = 0.0;

    for (int i = 1; i <= static_cast<int>(T / dividendtime); i++) {

        SumDiv += div * exp(-r * i * dividendtime);

    }

    S_adj = S0 - SumDiv;

    for (int i = 0; i < N + 1; i++) {

        double s = 0.0;

        for (int k = 1; k <= static_cast<int>(T / dividendtime); k++) {

            if (i * deltaT < k * dividendtime) {

                s += div * exp(-r * (k * dividendtime - i * deltaT));

            }

            for (int j = i; j >= 0; j--) {

                D[i][j] = 0.0 + s;

            }
        }
    }

    for (int j = N; j >= 0; j--) {

        S[N][j] = S_adj * pow(u, 2 * j - N) + D[N][j];
        V[N][j] = max(K - S[N][j], 0.0);

    }

    for (int i = N - 1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S_adj * pow(u, 2 * j - i) + D[i][j];
            V[i][j] = max(DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]), K - S[i][j]);

        }
    }

    double value =  V[0][0];
    double delta = (V[1][1] - V[1][0]) / (S[1][1] - S[1][0]);
    double gamma = ((V[2][2] - V[2][1]) / (S[2][2] - S[2][1]) - (V[2][1] - V[2][0]) / (S[2][1] - S[2][0])) / ((S[2][2] - S[2][0]) / 2);
    double theta = (V[2][1] - V[0][0]) / (2 * deltaT);

    return { value, delta, gamma, theta };

}


int main() {

    double spotPrice = 50.0;
    double strikePrice = 55.55;
    double volatility = 0.25;
    double dividendAmount = 0.95;
    double dividendPeriod = 2.0 / 12.0;
    double riskFreeRate = 0.04;
    double optionTerm = 7.0 / 12.0;

    int eurIniN = 7;
    int amerIniN = 7;

    double tolerance = pow(10, -4);

    double eurDiff = 1.0;
    double amerDiff = 1.0;

    vector <double> optEur_1;
    vector <double> optEur_2;
    vector <double> optAmer_1;
    vector <double> optAmer_2;

    
    while (eurDiff >= tolerance) {

            optEur_1 = EuropeanBinomialPutFixedDividend(eurIniN, spotPrice, strikePrice, riskFreeRate, dividendAmount, dividendPeriod, volatility, optionTerm);
            optEur_2 = EuropeanBinomialPutFixedDividend(eurIniN*2, spotPrice, strikePrice, riskFreeRate, dividendAmount, dividendPeriod, volatility, optionTerm);
            eurDiff = abs(optEur_1[0] - optEur_2[0]);
            eurIniN *= 2;
    }
    
    cout << optEur_1[0] << "\t" << eurIniN << "\t" << optEur_1[1] << "\t" << optEur_1[2] << "\t" << optEur_1[3] << endl;

    while (amerDiff >= tolerance) {

            optAmer_1 = AmericanBinomialPutFixedDividend(amerIniN, spotPrice, strikePrice, riskFreeRate, dividendAmount, dividendPeriod, volatility, optionTerm);
            optAmer_2 = AmericanBinomialPutFixedDividend(amerIniN * 2, spotPrice, strikePrice, riskFreeRate, dividendAmount, dividendPeriod, volatility, optionTerm);
            amerDiff = abs(optAmer_1[0] - optAmer_2[0]);
            amerIniN *= 2;
    }
     
    cout << optAmer_1[0] << "\t" << amerIniN << "\t" << optAmer_1[1] << "\t" << optAmer_1[2] << "\t" << optAmer_1[3] << endl;
    
    return 0;

}