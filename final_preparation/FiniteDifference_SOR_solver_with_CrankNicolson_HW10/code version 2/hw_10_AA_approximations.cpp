#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>


using namespace std;


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


vector <double> SOR_amer(vector <vector<double>> A, vector <double> B, vector <double> X0, vector <double> X_min, double w, double tol) {

	double norm = 10.0;
	int n = B.size();
	vector<double>X(n);

	while (norm >= tol) {

		double sumsquared = 0.0;

		for (int i = 0; i < n; i++) {

			double sum1 = 0.0;
			double sum2 = 0.0;

			for (int k = 0; k < i; k++) {

				sum1 += A[i][k] * X[k];

			}

			for (int k = i + 1; k < n; k++) {

				sum2 += A[i][k] * X0[k];

			}

			X[i] = max((1 - w) * X0[i] - w / A[i][i] * (sum1 + sum2) + w / A[i][i] * B[i], X_min[i]);
		}

		for (int i = 0; i < n; i++) {

			sumsquared += pow((X[i] - X0[i]), 2.0);

		}

		norm = sqrt(sumsquared);
		X0 = X;
	}

	return X;
}



vector <double> SOR_CrankNicolson_amer(double spot, double strike, double riskfree, double dividend, double sigma, double term, double m, double temp_alpha, vector <double> B, vector <double> X0, vector <double> X_min, double w, double tol) {

	double norm = 10.0;
	int n = B.size();
	vector<double>X(n);
	double x_left = log(spot / strike) + (riskfree - dividend - 0.5 * sigma * sigma) * term - 3.0 * sigma * sqrt(term);
	double x_right = log(spot / strike) + (riskfree - dividend - 0.5 * sigma * sigma) * term + 3.0 * sigma * sqrt(term);
	double tau_final = term * (sigma * sigma) / 2.0;
	double delta_tau = tau_final / m;
	double delta_t = 2.0 * delta_tau / (sigma * sigma);
	int nn = static_cast<int>((x_right - x_left) / sqrt(delta_tau / temp_alpha));
	double delta_x = (x_right - x_left) / nn;
	double alpha = delta_tau / (delta_x * delta_x);

	while (norm >= tol) {

		double sumsquared = 0.0;

		for (int i = 0; i < n; i++) {

			if (i == 0) {

				X[i] = max((1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * X0[i + 1] + w * B[i] / (1 + alpha), X_min[i]);

			}

			else {

				if (i == n - 1) {

					X[i] = max((1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * X[i - 1] + w * B[i] / (1 + alpha), X_min[i]);

				}

				else
				{

					X[i] = max((1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * (X[i - 1] + X0[i + 1]) + w * B[i] / (1 + alpha), X_min[i]);

				}
			}
		}

		for (int i = 0; i < n; i++) {

			sumsquared += pow((X[i] - X0[i]), 2.0);

		}

		norm = sqrt(sumsquared);
		X0 = X;

	}

	return X;

}


vector <vector<double>> Backward_CrankNicolson_SOR_FD_Put_American_Approximations(double spot, double strike, double riskfree, double dividend, double sigma, double term, double temp_alpha, double m, double w, double tol) {


	int beg_level = 0;

	double x_compute = log(spot / strike);
	double u_compute;

	double a = (riskfree - dividend) / (sigma * sigma) - 0.5;
	double b = pow((riskfree - dividend) / (sigma * sigma) + 0.5, 2.0) + 2.0 * dividend / (sigma * sigma);

	double tau_final = term * (sigma * sigma) / 2.0;
	double x_left = log(spot / strike) + (riskfree - dividend - 0.5 * sigma * sigma) * term - 3.0 * sigma * sqrt(term);
	double x_right = log(spot / strike) + (riskfree - dividend - 0.5 * sigma * sigma) * term + 3.0 * sigma * sqrt(term);

	double delta_tau = tau_final / m;
	double delta_t = 2.0 * delta_tau / (sigma * sigma);
	int n = static_cast<int>((x_right - x_left) / sqrt(delta_tau / temp_alpha));
	double delta_x = (x_right - x_left) / n;
	double alpha = delta_tau / (delta_x * delta_x);

	vector <double> x(n + 1);
	vector<double> tau(m + 1);
	vector <vector <double>> u_x(m + 1, vector<double>(n + 1));


	for (int i = 0; i <= m; i++) {

		tau[i] = i * delta_tau;

	}


	for (int j = 0; j <= n; j++) {

		x[j] = x_left + j * delta_x;

	}


	for (int i = 0; i <= m; i++) {

		u_x[i][0] = strike * exp(a * x_left + b * tau[i]) * (1 - exp(x_left));
		u_x[i][n] = 0.0;

	}


	for (int j = 1; j < n; j++) {

		u_x[0][j] = strike * exp(a * x[j]) * max(1 - exp(x[j]), 0.0);

	}


	vector <vector <double>> A(n - 1, vector<double>(n - 1));
	vector <vector <double>> AA(n - 1, vector<double>(n - 1));
	vector <double> B(n - 1);
	vector <double> X0(n - 1);
	vector <double> X_min(n - 1);
	vector <double> X(n - 1);

	for (int i = 0; i < n - 1; i++) {

		for (int j = 0; j < n - 1; j++) {

			if (i == j) {
				A[i][j] = 1 + alpha;
				AA[i][j] = 1 - alpha;
			}

			else {

				if (abs(i - j) == 1.0) {
					A[i][j] = -alpha / 2.0;
					AA[i][j] = alpha / 2.0;
				}

				else {
					A[i][j] = 0;
					AA[i][j] = 0;
				}
			}
		}
	}


	for (int i = 1; i <= m; i++) {

		vector <double> Temp(n - 1);
		for (int j = 0; j < n - 1; j++) {
			Temp[j] = u_x[i - 1][j + 1];
		}

		Temp = MatrixTimesVector(AA, Temp);

		B[0] = alpha / 2.0 * u_x[i][0] + alpha / 2.0 * u_x[i - 1][0] + Temp[0];
		B[n - 2] = alpha / 2.0 * u_x[i][n] + alpha / 2.0 * u_x[i - 1][n] + Temp[n - 2];

		for (int j = 1; j < n - 2; j++) {

			B[j] = Temp[j];

		}

		for (int j = 0; j < n - 1; j++) {

			X0[j] = strike * exp(a * x[j + 1] + b * tau[i]) * max(1 - exp(x[j + 1]), 0.0);

			X_min[j] = strike * exp(a * x[j + 1] + b * tau[i]) * max(1 - exp(x[j + 1]), 0.0);

		}

		//X = SOR_amer(A, B, X0, X_min, w, tol);

		X = SOR_CrankNicolson_amer(spot, strike, riskfree, dividend, sigma, term, m, temp_alpha, B, X0, X_min, w, tol);

		for (int j = 1; j < n; j++) {

			u_x[i][j] = X[j - 1];

		}

	}

	return u_x;

}
	


int main() {

	double spotPrice = 37.0;
	double strikePrice = 40.0;
	double riskFreeRate = 0.04;
	double dividendYield = 0.02;
	double Volatility = 0.28;
	double optionTerm = 0.75;
	char optionType = 'P';
	double P_amer_bin = 5.0455539623;

	double alpha_temp1 = 0.45;
	double alpha_temp2 = 5.0;
	int M = 4;

	double w = 1.2;
	double tol = pow(10, -6);

	vector <vector <double>> AmericanPut_FD_CN_SOR_Approximations = Backward_CrankNicolson_SOR_FD_Put_American_Approximations(spotPrice, strikePrice, riskFreeRate, dividendYield, Volatility, optionTerm, alpha_temp1, M, w, tol);


	for (int i = 0; i <= M; i++) {

		for (int j = 0; j < AmericanPut_FD_CN_SOR_Approximations[i].size(); j++) {
			cout << AmericanPut_FD_CN_SOR_Approximations[i][j] << "\t";
		}

		cout << endl;

	}
}

