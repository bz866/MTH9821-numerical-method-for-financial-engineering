# include <iostream>
# include <math.h>
# include <vector>

using namespace std;


vector <double> EuropeanBinomialPutDividend(int N, double S0, double K, double r, vector<double> fixedDividends, vector<double> fixedDividendsTerms, vector<double> percentageDividends, vector<double> percentageDividendsTerms, double sigma, double T) {

    vector <double> D_F(N + 1);
    vector <double> D_P(N + 1);
    vector<vector <double>> S(N + 1);
    vector<vector <double>> V(N + 1);


    for (int i = 0; i < N + 1 ; i++) {

        S[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    double sumDiv = 0.0;

    for (int i = 0; i < fixedDividends.size(); i++) {

        sumDiv += fixedDividends[i] * exp(-r * fixedDividendsTerms[i]);

    }
    
    double S_adj = S0 - sumDiv;


    for (int i = 0; i < N + 1; i++) {

        double s = 0.0;
        double prod = 1.0;

        for (int k = 0; k < fixedDividends.size(); k++) {

            if (i * deltaT < fixedDividendsTerms[k]) {

                s += fixedDividends[k] * exp(-r * (fixedDividendsTerms[k] - i * deltaT));
            }
        }


        for (int m = 0; m < percentageDividends.size(); m++) {

            if (i * deltaT >= percentageDividendsTerms[m]) {

                prod *= (1 - percentageDividends[m]); 
            }
        }

            D_F[i] = s;
            D_P[i] = prod;

    }
    

    for (int j = N; j >= 0; j--) {

        S[N][j] = S_adj * pow(u, 2 * j - N) * D_P[N] + D_F[N];
        V[N][j] = max(K - S[N][j], 0.0);

    }


    for (int i = N - 1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S_adj * pow(u, 2 * j - i) * D_P[i] + D_F[i];
            V[i][j] = DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]);

        }
    }

    double value = V[0][0];
    double delta = (V[1][1] - V[1][0]) / (S[1][1] - S[1][0]);
    double gamma = ((V[2][2] - V[2][1]) / (S[2][2] - S[2][1]) - (V[2][1] - V[2][0]) / (S[2][1] - S[2][0])) / ((S[2][2] - S[2][0]) / 2);
    double theta = (V[2][1] - V[0][0]) / (2 * deltaT);

    return { value, delta, gamma, theta };

}


vector <double> AmericanBinomialPutDividend(int N, double S0, double K, double r, vector<double> fixedDividends, vector<double> fixedDividendsTerms, vector<double> percentageDividends, vector<double> percentageDividendsTerms, double sigma, double T) {

    vector <double> D_F(N + 1);
    vector <double> D_P(N + 1);
    vector<vector <double>> S(N + 1);
    vector<vector <double>> V(N + 1);


    for (int i = 0; i < N + 1; i++) {

        S[i].resize(i + 1);
        V[i].resize(i + 1);

    }

    double deltaT = T / N;
    double u = exp(sigma * sqrt(deltaT));
    double d = 1 / u;

    double P_u = (exp(r * deltaT) - d) / (u - d);
    double P_d = 1 - P_u;
    double DF = exp(-r * deltaT);

    double sumDiv = 0.0;

    for (int i = 0; i < fixedDividends.size(); i++) {

        sumDiv += fixedDividends[i] * exp(-r * fixedDividendsTerms[i]);

    }

    double S_adj = S0 - sumDiv;


    for (int i = 0; i < N + 1; i++) {

        double s = 0.0;
        double prod = 1.0;

        for (int k = 0; k < fixedDividends.size(); k++) {

            if (i * deltaT < fixedDividendsTerms[k]) {

                s += fixedDividends[k] * exp(-r * (fixedDividendsTerms[k] - i * deltaT));
            }
        }


        for (int m = 0; m < percentageDividends.size(); m++) {

            if (i * deltaT >= percentageDividendsTerms[m]) {

                prod *= (1 - percentageDividends[m]);
            }
        }

        D_F[i] = s;
        D_P[i] = prod;

    }


    for (int j = N; j >= 0; j--) {

        S[N][j] = S_adj * pow(u, 2 * j - N) * D_P[N] + D_F[N];
        V[N][j] = max(K - S[N][j], 0.0);

    }


    for (int i = N - 1; i >= 0; i--) {

        for (int j = i; j >= 0; j--) {

            S[i][j] = S_adj * pow(u, 2 * j - i) * D_P[i] + D_F[i];
            V[i][j] = max(DF * (P_u * V[i + 1][j + 1] + P_d * V[i + 1][j]),K- S[i][j]);

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
    double riskFreeRate = 0.04;
    double optionTerm = 7.0 / 12.0;

    vector<double> fixedDividends{0.95, 0.45};
    vector<double> fixedDividendsTerms{2.0/12.0, 6.0/12.0};

    vector<double> percentageDividends{0.02};
    vector<double> percentageDividendsTerms{4.0/12.0};

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

        optEur_1 = EuropeanBinomialPutDividend(eurIniN, spotPrice, strikePrice, riskFreeRate, fixedDividends, fixedDividendsTerms, percentageDividends, percentageDividendsTerms, volatility, optionTerm);
        optEur_2 = EuropeanBinomialPutDividend(eurIniN*2, spotPrice, strikePrice, riskFreeRate, fixedDividends, fixedDividendsTerms, percentageDividends, percentageDividendsTerms, volatility, optionTerm);
        eurDiff = abs(optEur_1[0] - optEur_2[0]);
        eurIniN *= 2;
    }

    cout << optEur_1[0] << "\t" << eurIniN << "\t" << optEur_1[1] << "\t" << optEur_1[2] << "\t" << optEur_1[3] << endl;

    while (amerDiff >= tolerance) {

        optAmer_1 = AmericanBinomialPutDividend(amerIniN, spotPrice, strikePrice, riskFreeRate, fixedDividends, fixedDividendsTerms, percentageDividends, percentageDividendsTerms, volatility, optionTerm);
        optAmer_2 = AmericanBinomialPutDividend(amerIniN * 2, spotPrice, strikePrice, riskFreeRate, fixedDividends, fixedDividendsTerms, percentageDividends, percentageDividendsTerms, volatility, optionTerm);
        amerDiff = abs(optAmer_1[0] - optAmer_2[0]);
        amerIniN *= 2;
    }

    cout << optAmer_1[0] << "\t" << amerIniN << "\t" << optAmer_1[1] << "\t" << optAmer_1[2] << "\t" << optAmer_1[3] << endl;
    
    
    return 0;

}