#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>

using namespace std;


vector <double> ForwardEulerFD_Call_European_PercentageDividend_Greeks(double spot, double strike, double riskfree, double dividendPercentage, double dividendTime, double sigma, double term, double alpha1, double M1) {


	double a = riskfree / (sigma * sigma) - 0.5;
	double b = pow(riskfree / (sigma * sigma) + 0.5, 2.0);

	double tau_final = term * (sigma * sigma) / 2.0;
	double tau_div = (term - dividendTime) * (sigma * sigma) / 2.0;
	double delta_tau_1 = tau_div / M1;
	double delta_x = sqrt(delta_tau_1 / alpha1);

	double x_compute_bar = log(spot / strike) + log(1 - dividendPercentage);
	double x_temp_left = log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term - 3.0 * sigma * sqrt(term);
	double x_temp_right = log(spot / strike) + (riskfree - 0.5 * sigma * sigma) * term + 3.0 * sigma * sqrt(term);

	double N_left = ceil((x_compute_bar - x_temp_left) / delta_x);
	double N_right = ceil((x_temp_right - x_compute_bar) / delta_x);
	double N = N_left + N_right;

	double x_left = x_compute_bar - N_left * delta_x;
	double x_right = x_compute_bar + N_right * delta_x;

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
	double M2 = ceil((tau_final - tau_div) / delta_tau_2_temp);
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

	double x_compute = x_compute_bar - log(1 - dividendPercentage);
	double x_low = x_compute - delta_x;
	double x_high = x_compute + delta_x;

	double s = strike * exp(x_compute);
	double s_low = strike * exp(x_low);
	double s_high = strike * exp(x_high);

	double value_0 = exp(-a * x_compute - b * tau_final) * u2_x[M2][N_left];
	double value_0_low = exp(-a * x_low - b * tau_final) * u2_x[M2][N_left - 1];
	double value_0_high = exp(-a * x_high - b * tau_final) * u2_x[M2][N_left + 1];
	double value_1 = exp(-a * x_compute - b * (tau_final - delta_tau_2)) * u2_x[M2 - 1][N_left];
	
	double delta = (value_0_high - value_0_low) / (s_high - s_low);
	double gamma = ((value_0_high - value_0) / (s_high - s) - (value_0 - value_0_low) / (s - s_low)) / ((s_high - s_low) / 2);
	double theta = (value_1 - value_0) / (delta_tau_2 / (0.5 * sigma * sigma));

	return{M2, alpha2, N, x_left, x_right, x_left_new, x_right_new, tau_div, delta_tau_1, delta_tau_2, delta_x, u2_x[M2][N_left], value_0, delta, gamma, theta};

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

	for (int i = 1; i <= 4; i++) {

		vector<double> FD_Call = ForwardEulerFD_Call_European_PercentageDividend_Greeks(spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendTime, Volatility, optionTerm, alpha1, pow(4.0, i));

		cout << FD_Call[0] << "\t" << FD_Call[1] << "\t" << FD_Call[2] << "\t" << FD_Call[3] << "\t" << FD_Call[4] << "\t" << FD_Call[5] << "\t" << FD_Call[6] << "\t" << FD_Call[7] << "\t" << FD_Call[8]
			<< "\t" << FD_Call[9] << "\t" << FD_Call[10] << "\t" << FD_Call[11] << "\t" << FD_Call[12] << "\t" << FD_Call[13] << "\t" << FD_Call[14] << "\t" << FD_Call[15] << endl;
	}


}