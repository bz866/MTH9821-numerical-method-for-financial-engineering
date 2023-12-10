#include<iostream>
#include<math.h>
#include<iomanip>
#include<vector>


using namespace std;

double BSM(double S0, double K, double r, double y, double sigma, double T, char type) {

	double value;

	double Nd1;
	double Nd2;

	double d1 = (log(S0 / K) + (r - y + .5 * sigma * sigma) * T) / (sigma * sqrt(T));
	double d2 = d1 - sigma * sqrt(T);

	Nd1 = .5 * (1 + erf(d1 / sqrt(2.0)));
	Nd2 = .5 * (1 + erf(d2 / sqrt(2.0)));

	value = Nd1 * S0 * exp(-y * T) - Nd2 * K * exp(-r * T);

	if (type == 'P') {

		value = value - S0 * exp(-y * T) + K * exp(-r * T);

	}

	return value;

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

			sumsquared += pow(X[i] - X0[i], 2.0);

		}

		norm = sqrt(sumsquared);
		X0 = X;
	}

	return X;

}

vector <double> SOR_eur(vector <vector<double>> A, vector <double> B, vector <double> X0, double w, double tol) {

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

			X[i] = (1 - w) * X0[i] - w / A[i][i] * (sum1 + sum2) + w / A[i][i] * B[i];
		}

		for (int i = 0; i < n; i++) {

			sumsquared += pow(X[i] - X0[i], 2.0);

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

			sumsquared += pow(X[i] - X0[i], 2.0);

		}

		norm = sqrt(sumsquared);
		X0 = X;

	}

	return X;

}

vector <double> SOR_CrankNicolson_eur(double spot, double strike, double riskfree, double dividend, double sigma, double term, double m, double temp_alpha, vector <double> B, vector <double> X0, double w, double tol) {

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

				X[i] = (1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * X0[i + 1] + w * B[i] / (1 + alpha);

			}

			else {

				if (i == n - 1) {

					X[i] = (1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * X[i - 1] + w * B[i] / (1 + alpha);

				}

				else
				{

					X[i] = (1 - w) * X0[i] + w * alpha / (2 * (1 + alpha)) * (X[i - 1] + X0[i + 1]) + w * B[i] / (1 + alpha);

				}
			}
		}

		for (int i = 0; i < n; i++) {

			sumsquared += pow(X[i] - X0[i], 2.0);

		}

		norm = sqrt(sumsquared);
		X0 = X;

	}

	return X;

}


vector<double> Backward_CrankNicolson_SOR_FD_Put_European_Value(double spot, double strike, double riskfree, double dividend, double sigma, double term, double temp_alpha, double m, double w, double tol) {

	int beg_level = 0;

	double BSM_put = BSM(spot, strike, riskfree, dividend, sigma, term, 'P');
	double x_compute = log(spot / strike);

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

	double value_approx_eur;

	double x_low;
	double x_high;

	double s_low;
	double s_high;

	double value_0_low;
	double value_0_high;

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

		u_x[i][0] = strike * exp(a * x_left + b * tau[i]) * (exp(-2 * riskfree * tau[i] / (sigma * sigma)) - exp(x_left - 2 * dividend * tau[i] / (sigma * sigma)));
		u_x[i][n] = 0.0;

	}


	for (int j = 1; j < n; j++) {

		u_x[0][j] = strike * exp(a * x[j]) * max(1 - exp(x[j]), 0.0);

	}


	vector <vector <double>> A(n - 1, vector<double>(n - 1));
	vector <vector <double>> AA(n - 1, vector<double>(n - 1));
	vector <double> B(n - 1);
	vector <double> X0(n - 1);
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

		}


		X = SOR_CrankNicolson_eur(spot, strike, riskfree, dividend, sigma, term, m, temp_alpha, B, X0, w, tol);


		for (int j = 1; j < n; j++) {

			u_x[i][j] = X[j - 1];
		}

	}


	while (x[beg_level] <= x_compute) {
		beg_level++;
	}
	beg_level--;

	x_low = x[beg_level];
	x_high = x[beg_level + 1];

	s_low = strike * exp(x_low);
	s_high = strike * exp(x_high);

	value_0_low = exp(-a * x_low - b * tau_final) * u_x[m][beg_level];
	value_0_high = exp(-a * x_high - b * tau_final) * u_x[m][beg_level + 1];

	value_approx_eur = ((s_high - spot) * value_0_low + (spot - s_low) * value_0_high) / (s_high - s_low);

	return {BSM_put, value_approx_eur};

}



	vector<double> Backward_CrankNicolson_SOR_FD_Put_American_Greeks_Errors(double spot, double strike, double riskfree, double dividend, double sigma, double term, double exactVal, double temp_alpha, double m, double w, double tol) {

		int beg_level = 0;

		double delta;
		double gamma;
		double theta;

		double error_pointwise_1;
		double error_pointwise_2;

		double BSM_put = BSM(spot, strike, riskfree, dividend, sigma, term, 'P');
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

		double value_approx_1;
		double value_approx_2;
		double value_approx_next;

		double x_low_low;
		double x_low;
		double x_high;
		double x_high_high;
		double s_low_low;
		double s_low;
		double s_high;
		double s_high_high;
		double value_0_low_low;
		double value_0_low;
		double value_0_high;
		double value_0_high_high;
		double value_1_low;
		double value_1_high;

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
			B[n-2] = alpha / 2.0 * u_x[i][n] + alpha / 2.0 * u_x[i - 1][n] + Temp[n-2];

			for (int j = 1; j < n - 2; j++) {

			    B[j] = Temp[j];
							
			}

			for (int j = 0; j < n - 1; j++) {

				X0[j] = strike * exp(a * x[j + 1] + b * tau[i]) * max(1 - exp(x[j + 1]), 0.0);

				X_min[j] = strike * exp(a * x[j + 1] + b * tau[i]) * max(1 - exp(x[j + 1]), 0.0);
			}

			X = SOR_CrankNicolson_amer(spot, strike, riskfree, dividend, sigma, term, m, temp_alpha, B, X0, X_min, w, tol);

			for (int j = 1; j < n; j++) {

				u_x[i][j] = X[j - 1];

			}

		}


		while (x[beg_level] <= x_compute) {
			beg_level++;
		}
		beg_level--;

		x_low_low = x[beg_level - 1];
		x_low = x[beg_level];
		x_high = x[beg_level + 1];
		x_high_high = x[beg_level + 2];

		s_low_low = strike * exp(x_low_low);
		s_low = strike * exp(x_low);
		s_high = strike * exp(x_high);
		s_high_high = strike * exp(x_high_high);

		value_0_low_low = exp(-a * x_low_low - b * tau_final) * u_x[m][beg_level - 1];
		value_0_low = exp(-a * x_low - b * tau_final) * u_x[m][beg_level];
		value_0_high = exp(-a * x_high - b * tau_final) * u_x[m][beg_level + 1];
		value_0_high_high = exp(-a * x_high_high - b * tau_final) * u_x[m][beg_level + 2];

		value_approx_1 = ((s_high - spot) * value_0_low + (spot - s_low) * value_0_high) / (s_high - s_low);
		error_pointwise_1 = abs(value_approx_1 - exactVal);

		u_compute = ((x_high - x_compute) * u_x[m][beg_level] + (x_compute - x_low) * u_x[m][beg_level + 1]) / (x_high - x_low);
		value_approx_2 = exp(-a * x_compute - b * tau_final) * u_compute;
		error_pointwise_2 = abs(value_approx_2 - exactVal);

		delta = (value_0_high - value_0_low) / (s_high - s_low);
		gamma = ((value_0_high_high - value_0_high) / (s_high_high - s_high) - (value_0_low - value_0_low_low) / (s_low - s_low_low)) / ((s_high_high + s_high) / 2 - (s_low + s_low_low) / 2);

		value_1_low = exp(-a * x_low - b * (tau_final - delta_tau)) * u_x[m - 1][beg_level];
		value_1_high = exp(-a * x_high - b * (tau_final - delta_tau)) * u_x[m - 1][beg_level + 1];
		value_approx_next = ((s_high - spot) * value_1_low + (spot - s_low) * value_1_high) / (s_high - s_low);

		theta = (value_approx_next - value_approx_1) / delta_t;


		return { error_pointwise_1 ,error_pointwise_2, delta, gamma, theta, value_approx_1, exactVal};

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

	cout << "Alpha =" << alpha_temp1 << endl;

	for (int i = 1; i <= 4; i++) {

		vector<double> FD_Put_Amer = Backward_CrankNicolson_SOR_FD_Put_American_Greeks_Errors(spotPrice, strikePrice, riskFreeRate, dividendYield, Volatility, optionTerm, P_amer_bin, alpha_temp1, pow(4,i), w, tol);
		vector<double> FD_Put_Eur = Backward_CrankNicolson_SOR_FD_Put_European_Value(spotPrice, strikePrice, riskFreeRate, dividendYield, Volatility, optionTerm, alpha_temp1, pow(4, i), w, tol);

		cout << FD_Put_Amer[0] << "\t" << FD_Put_Amer[1] << "\t" << FD_Put_Amer[2] << "\t" << FD_Put_Amer[3] << "\t" << FD_Put_Amer[4] << "\t" << FD_Put_Amer[5]+ FD_Put_Eur[0]-FD_Put_Eur[1] << "\t" << abs(FD_Put_Amer[5] + FD_Put_Eur[0] - FD_Put_Eur[1]- FD_Put_Amer[6]) << endl;

	}

	cout << endl;
	cout << endl;
	
	cout << "Alpha =" << alpha_temp2 << endl;
	for (int i = 1; i <= 4; i++) {

		vector<double> FD_Put_Amer = Backward_CrankNicolson_SOR_FD_Put_American_Greeks_Errors(spotPrice, strikePrice, riskFreeRate, dividendYield, Volatility, optionTerm, P_amer_bin, alpha_temp2, pow(4,i), w, tol);
		vector<double> FD_Put_Eur = Backward_CrankNicolson_SOR_FD_Put_European_Value(spotPrice, strikePrice, riskFreeRate, dividendYield, Volatility, optionTerm, alpha_temp2, pow(4, i), w, tol);

		cout << FD_Put_Amer[0] << "\t" << FD_Put_Amer[1] << "\t" << FD_Put_Amer[2] << "\t" << FD_Put_Amer[3] << "\t" << FD_Put_Amer[4] << "\t" << FD_Put_Amer[5] + FD_Put_Eur[0] - FD_Put_Eur[1] << "\t" << abs(FD_Put_Amer[5] + FD_Put_Eur[0] - FD_Put_Eur[1] - FD_Put_Amer[6]) << endl;


	}

}