#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>


using namespace std;


vector <vector <double>> ForwardEulerFD_Call_European_PercentageDividend_Approximations1(double spot, double strike, double riskfree, double dividendPercentage, double dividendTime, double sigma, double term, double alpha1, int M1) {


	double a = riskfree / (sigma * sigma) - 0.5;
	double b = pow(riskfree / (sigma * sigma) + 0.5, 2.0);

	double tau_final = term * (sigma * sigma) / 2.0;
	double tau_div = (term - dividendTime) * (sigma * sigma) / 2.0;
	double delta_tau_1 = tau_div / M1;
	double delta_x = sqrt(delta_tau_1 / alpha1);

	double x_compute = log(spot / strike) + log(1 - dividendPercentage);
	double x_temp_left= log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term - 3.0 * sigma * sqrt(term);
	double x_temp_right = log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term + 3.0 * sigma * sqrt(term);

	int N_left = ceil((x_compute - x_temp_left) / delta_x);
	int N_right = ceil((x_temp_right - x_compute) / delta_x);
	int N = N_left + N_right;

	double x_left = x_compute - N_left * delta_x;
	double x_right = x_compute + N_right * delta_x;

	vector <double> x(N + 1);
	vector<double> tau_1(M1 + 1);
	vector <vector <double>> u1_x(M1 + 1, vector<double>(N + 1));


	for (int i = 0; i <= M1; i++) {
		tau_1[i] = i * delta_tau_1;
	}


	for (int j = 0; j <= N; j++) {
		x[j] = x_left + j * delta_x;
	}


	for (int i = 0; i <= M1; i++) {

		u1_x[i][0] = 0;
		u1_x[i][N] = strike * exp(a * x_right + b * tau_1[i]) * (exp(x_right) - exp((-2 * riskfree * tau_1[i]) / (sigma * sigma)));
	}


	for (int j = 1; j < N; j++) {

		u1_x[0][j] = strike * exp(a * x[j]) * max(exp(x[j]) - 1, 0.0);

	}


	for (int i = 1; i <= M1; i++) {

		for (int j = 1; j < N; j++) {

			u1_x[i][j] = alpha1 * u1_x[i - 1][j - 1] + (1 - 2 * alpha1) * u1_x[i - 1][j] + alpha1 * u1_x[i - 1][j + 1];

		}

	}

	return u1_x;

}



vector <vector <double>> ForwardEulerFD_Call_European_PercentageDividend_Approximations2(double spot, double strike, double riskfree, double dividendPercentage, double dividendTime, double sigma, double term, double alpha1, double M1) {


	double a = riskfree / (sigma * sigma) - 0.5;
	double b = pow(riskfree / (sigma * sigma) + 0.5, 2.0);

	double tau_final = term * (sigma * sigma) / 2.0;
	double tau_div = (term - dividendTime) * (sigma * sigma) / 2.0;
	double delta_tau_1 = tau_div / M1;
	double delta_x = sqrt(delta_tau_1 / alpha1);

	double x_compute = log(spot / strike) + log(1 - dividendPercentage);
	double x_temp_left = log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term - 3.0 * sigma * sqrt(term);
	double x_temp_right = log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term + 3.0 * sigma * sqrt(term);

	int N_left = ceil((x_compute - x_temp_left) / delta_x);
	int N_right = ceil((x_temp_right - x_compute) / delta_x);
	int N = N_left + N_right;

	double x_left = x_compute - N_left * delta_x;
	double x_right = x_compute + N_right * delta_x;

	vector <double> x(N + 1);
	vector<double> tau_1(M1 + 1);
	vector <vector <double>> u1_x(M1 + 1, vector<double>(N + 1));


	for (int i = 0; i <= M1; i++) {
		tau_1[i] = i * delta_tau_1;
	}


	for (int j = 0; j <= N; j++) {
		x[j] = x_left + j * delta_x;
	}


	for (int i = 0; i <= M1; i++) {

		u1_x[i][0] = 0;
		u1_x[i][N] = strike * exp(a * x_right + b * tau_1[i]) * (exp(x_right) - exp((-2 * riskfree * tau_1[i]) / (sigma * sigma)));
	}


	for (int j = 1; j < N; j++) {

		u1_x[0][j] = strike * exp(a * x[j]) * max(exp(x[j]) - 1, 0.0);

	}


	for (int i = 1; i <= M1; i++) {

		for (int j = 1; j < N; j++) {

			u1_x[i][j] = alpha1 * u1_x[i - 1][j - 1] + (1 - 2 * alpha1) * u1_x[i - 1][j] + alpha1 * u1_x[i - 1][j + 1];

		}

	}

	double alpha_temp = alpha1;
	double delta_tau_2_temp = alpha_temp * delta_x * delta_x;
	int M2 = ceil((tau_final - tau_div) / delta_tau_2_temp);
	double delta_tau_2 = (tau_final - tau_div) / M2;
	double alpha2 = delta_tau_2 / (delta_x * delta_x);
	double x_left_new = x_left - log(1 - dividendPercentage);
	double x_right_new = x_right - log(1 - dividendPercentage);


	vector <double> x_new(N + 1);
	vector<double> tau_2(M2 + 1);
	vector <vector <double>> u2_x(M2 + 1, vector<double>(N + 1));

	for (int i = 0; i <= M2; i++) {
		tau_2[i] = tau_div + i * delta_tau_2;
	}

	for (int j = 0; j <= N; j++) {
		x_new[j] = x_left_new + j * delta_x;
	}

	for (int i = 0; i <= M2; i++) {

		u2_x[i][0] = 0;
		u2_x[i][N] = strike * exp(a * x_right_new + b * tau_2[i]) * (exp(x_right) - exp((-2 * riskfree * tau_2[i]) / (sigma * sigma)));
	}

	for (int j = 1; j < N; j++) {

		u2_x[0][j] = pow(1 - dividendPercentage, -a) * u1_x[M1][j];

	}


	for (int i = 1; i <= M2; i++) {

		for (int j = 1; j < N; j++) {

			u2_x[i][j] = alpha2 * u2_x[i - 1][j - 1] + (1 - 2 * alpha2) * u2_x[i - 1][j] + alpha2 * u2_x[i - 1][j + 1];

		}

	}

	return u2_x;

}



int main() {

	double spotPrice = 52.0;
	double strikePrice = 50.0;
	double riskFreeRate = 0.03;
	double dividendPercentage = 0.01;
	double dividendTime = 5.0 / 12.0;
	double Volatility = 0.20;
	double optionTerm = 1.0;
	char optionType = 'C';

	double alpha1 = 0.4;
	int M1 = 4;

	vector <vector <double>> EuropeanPutDividends_FD_Approximations1 = ForwardEulerFD_Call_European_PercentageDividend_Approximations1(spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendTime, Volatility, optionTerm, alpha1, M1);

	for (int i = 0; i <= M1; i++) {

		for (int j = 0; j < EuropeanPutDividends_FD_Approximations1[i].size(); j++) {
			cout << EuropeanPutDividends_FD_Approximations1[i][j] << "\t";
		}

		cout << endl;

	}

	cout << endl;

	vector <vector <double>> EuropeanPutDividends_FD_Approximations2 = ForwardEulerFD_Call_European_PercentageDividend_Approximations2(spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendTime, Volatility, optionTerm, alpha1, M1);
	
		for (int i = 0; i < EuropeanPutDividends_FD_Approximations2.size(); i++) {

			for (int j = 0; j < EuropeanPutDividends_FD_Approximations2[i].size(); j++) {
				cout << EuropeanPutDividends_FD_Approximations2[i][j] << "\t";
			}

			cout << endl;

		}
     
		
}