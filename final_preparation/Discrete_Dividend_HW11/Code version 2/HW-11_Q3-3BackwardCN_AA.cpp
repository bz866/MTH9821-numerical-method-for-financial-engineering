#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>

using namespace std;


vector <double> forwardL(vector <vector<double>> L, vector <double> B) {

	int n = static_cast<int>(B.size());
	vector <double> y(n);

	y[0] = B[0];
	for (int i = 1; i < n; i++) {
		y[i] = B[i] - L[i][i - 1] * y[i - 1];
	}

	return y;
}

vector <double> BackwardU(vector <vector<double>> U, vector <double> Y) {

	int n = static_cast<int>(Y.size());
	vector <double> x(n);

	x[n - 1] = Y[n - 1] / U[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		x[i] = (Y[i] - U[i][i + 1] * x[i + 1]) / U[i][i];
	}

	return x;
}


vector <double> tridiag_matrix_bidiag_LU(vector <vector<double>> A, vector <double> B) {

	int n = static_cast<int>(B.size());
	vector <double> x(n);
	vector <double> y(n);
	vector <vector <double>> L(n, vector<double>(n));
	vector <vector <double>> U(n, vector<double>(n));

	for (int i = 0; i < n - 1; i++) {
		L[i][i] = 1;
		L[i + 1][i] = A[i + 1][i] / A[i][i];
		U[i][i] = A[i][i];
		U[i][i + 1] = A[i][i + 1];
		A[i + 1][i + 1] = A[i + 1][i + 1] - L[i + 1][i] * U[i][i + 1];
	}

	L[n - 1][n - 1] = 1;
	U[n - 1][n - 1] = A[n - 1][n - 1];

	y = forwardL(L, B);
	x = BackwardU(U, y);

	return x;
}

vector <double> MatrixTimesVector(vector <vector<double>> A, vector <double> B) {

	int n = static_cast<int>(B.size());
	vector <double> x(n);

	for (int i = 0; i < n; i++) {
		double sum = 0;
		for (int j = 0; j < n; j++) {
			sum += A[i][j] * B[j];
		}
		x[i] = sum;
	}

	return x;
}

vector <double> BackwardCrankNicolsonFD_Call_European_PercentageDividend_Greeks(double spot, double strike, double riskfree, double dividendPercentage, double dividendTime, double sigma, double term, double alpha1, double M1) {


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

	
    vector <vector <double>> A(N - 1, vector<double>(N - 1));
	vector <vector <double>> AA(N - 1, vector<double>(N - 1));
	vector <double> B(N - 1);
	vector <double> X(N - 1);

	for (int i = 0; i < N - 1; i++) {

		for (int j = 0; j < N - 1; j++) {
			if (i == j) {
				A[i][j] = 1 + alpha1;
				AA[i][j] = 1 - alpha1;
			}

			else {

				if (abs(i - j) == 1.0) {
					A[i][j] = -alpha1 / 2.0;
					AA[i][j] = alpha1 / 2.0;
				}

				else {
					A[i][j] = 0;
					AA[i][j] = 0;
				}
			}
		}
	}



	for (int i = 1; i <= M1; i++) {

		vector <double> Temp(N - 1);
		for (int j = 0; j < N - 1; j++) {
			Temp[j] = u1_x[i - 1][j + 1];
		}

		Temp = MatrixTimesVector(AA, Temp);

		for (int j = 0; j < N - 1; j++) {

			if (j == 0) {
				B[j] = alpha1 / 2.0 * u1_x[i][0] + alpha1 / 2.0 * u1_x[i - 1][0] + Temp[j];
			}

			else {

				if (j == N - 2) {
					B[j] = alpha1 / 2.0 * u1_x[i][N] + alpha1 / 2.0 * u1_x[i - 1][N] + Temp[j];
				}

				else {
					B[j] = Temp[j];
				}
			}
		}

		X = tridiag_matrix_bidiag_LU(A, B);

		for (int j = 1; j < N; j++) {

			u1_x[i][j] = X[j - 1];

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


	
	vector <vector <double>> A_new(N - 1, vector<double>(N - 1));
	vector <vector <double>> AA_new(N - 1, vector<double>(N - 1));
	vector <double> B_new(N - 1);
	vector <double> X_new(N - 1);

	for (int i = 0; i < N - 1; i++) {

		for (int j = 0; j < N - 1; j++) {
			if (i == j) {
				A_new[i][j] = 1 + alpha2;
				AA_new[i][j] = 1 - alpha2;
			}

			else {

				if (abs(i - j) == 1.0) {
					A_new[i][j] = -alpha2 / 2.0;
					AA_new[i][j] = alpha2 / 2.0;
				}

				else {
					A_new[i][j] = 0;
					AA_new[i][j] = 0;
				}
			}
		}
	}



	for (int i = 1; i <= M2; i++) {

		vector <double> Temp_new(N - 1);
		for (int j = 0; j < N - 1; j++) {
			Temp_new[j] = u2_x[i - 1][j + 1];
		}

		Temp_new = MatrixTimesVector(AA_new, Temp_new);

		for (int j = 0; j < N - 1; j++) {

			if (j == 0) {
				B_new[j] = alpha2 / 2.0 * u2_x[i][0] + alpha2 / 2.0 * u2_x[i - 1][0] + Temp_new[j];
			}

			else {

				if (j == N - 2) {
					B_new[j] = alpha2 / 2.0 * u2_x[i][N] + alpha2 / 2.0 * u2_x[i - 1][N] + Temp_new[j];
				}

				else {
					B_new[j] = Temp_new[j];
				}
			}
		}

		X_new = tridiag_matrix_bidiag_LU(A_new, B_new);

		for (int j = 1; j < N; j++) {

			u2_x[i][j] = X_new[j - 1];

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

	return{ M2, alpha2, N, x_left, x_right, x_left_new, x_right_new, tau_div, delta_tau_1, delta_tau_2, delta_x, u2_x[M2][N_left], value_0, delta, gamma, theta };

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

	double alpha1 = 4.0;

	for (int i = 1; i <= 4; i++) {

		vector<double> FD_Call = BackwardCrankNicolsonFD_Call_European_PercentageDividend_Greeks(spotPrice, strikePrice, riskFreeRate, dividendPercentage, dividendTime, Volatility, optionTerm, alpha1, pow(4.0, i));

		cout << FD_Call[0] << "\t" << FD_Call[1] << "\t" << FD_Call[2] << "\t" << FD_Call[3] << "\t" << FD_Call[4] << "\t" << FD_Call[5] << "\t" << FD_Call[6] << "\t" << FD_Call[7] << "\t" << FD_Call[8]
			<< "\t" << FD_Call[9] << "\t" << FD_Call[10] << "\t" << FD_Call[11] << "\t" << FD_Call[12] << "\t" << FD_Call[13] << "\t" << FD_Call[14] << "\t" << FD_Call[15] << endl;
	}


}

