#include"Tree.hpp"
#include"BS.hpp"
#include<cmath>

// return max of two numbers 
double max(double d1, double d2)
{
	if (d1 > d2)
	{
		return d1; 
	}
	else
	{
		return d2;
	}
}

// using Binomial Tree to price options 
vector<double> BT_euro(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1.0 / u;
	double P_up = (exp((r - q) * dt) - d) / (u - d);
	double P_down = 1 - P_up;
	
	vector<double> nodes;
	double d_square = d * d;
	double S_temp = S * pow(u, N);
	
	for (int i = 0; i < N + 1; i++)
	{
		double S_T = S_temp;
		S_temp *= d_square;
		double V_T = 0;
		if (option == 'C')
		{
			V_T = max(S_T - K, 0);
			
		}
		else
		{
			V_T = max(K-S_T, 0);
		}
		nodes.push_back(V_T);
	}

	for (int i = N - 1; i > 1; i--)
	{
		for (int j = 0; j < i + 1; j++)
		{
			nodes[j] = (P_up * nodes[j] + P_down * nodes[j + 1]);
		}
	}

	double S_10 = S * u;
	double S_11 = S * d;

	double S_20 = S * u * u;
	double S_21 = S;
	double S_22 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];

	double V_10 = (P_up * V_20 + P_down * V_21);
	double V_11 =(P_up * V_21 + P_down * V_22);

	double V = exp(-r * T) * (P_up * V_10 + P_down * V_11);

	double delta = exp(-r * T / N * (N - 1)) * (V_10 - V_11) / (S_10 - S_11);
	double delta_10 = (V_20 - V_21) / (S_20 - S_21);
	double delta_11 = (V_21 - V_22) / (S_21 - S_22);
	dt = (double)T / N;
	double theta = (V_21* exp(-r * T/ (double)N* (double)(N-2)) - V) / (2 * dt);

	double gamma = 2 * exp(-r * T / N * (N - 1)) * (delta_10 - delta_11) / (S_20 - S_22);

	return vector<double>{V, delta, theta, gamma};
}


// using Binomial Tree to price american options 
vector<double> BT_am(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1.0 / u;
	double P_up = (exp((r - q) * dt) - d) / (u - d);
	double P_down = 1 - P_up;
	

	vector<double> nodes;
	vector<double> S_list;
	double d_square = d * d;
	double S_temp = S * pow(u, N);
	double disc = exp(-r * dt);

	for (int i = 0; i < N + 1; i++)
	{
		double S_T = S_temp;
		
		double V_T = 0;
		if (option == 'C')
		{
			V_T = max(S_T - K, 0);
		}
		else
		{
			V_T = max(K - S_T, 0);
		}
		nodes.push_back(V_T);
		S_list.push_back(S_temp);
		S_temp *= d_square;
	}

	for (int i = N - 1; i > 1; i--)
	{
		for (int j = 0; j < i + 1; j++)
		{
			S_list[j] = S_list[j] * d;
			auto node_temp1 = disc * (P_up * nodes[j] + P_down * nodes[j + 1]);
			double node_temp2 = 0;
			if (option == 'C')
			{
				node_temp2 = max(S_list[j] - K, 0);
			}
			else
			{
				node_temp2 = max(K - S_list[j], 0);
			}
			
			nodes[j] = max(node_temp1, node_temp2);
		}
	
	}

	double S_10 = S * u;
	double S_11 = S * d;

	double S_20 = S * u * u;
	double S_21 = S;
	double S_22 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];

	double V_10am = 0;
	double V_11am = 0;
	if (option == 'C')
	{
		V_10am = max(S_10 - K, 0);
		V_11am = max(S_11 - K, 0);
	}
	else
	{
		V_10am = max(K - S_10, 0);
		V_11am = max(K - S_11, 0);
	}
	double V_10 = max(disc * (P_up * V_20 + P_down * V_21), V_10am);
	double V_11 = max(disc * (P_up * V_21 + P_down * V_22), V_10am);

	double V_am = 0;
	if (option == 'C')
	{
		V_am = max(S - K, 0);
	}
	else
	{
		V_am = max(K - S, 0);
	}
	double V = max(disc * (P_up * V_10 + P_down * V_11), V_am);

	double delta = (V_10 - V_11) / (S_10 - S_11);
	double delta_10 = (V_20 - V_21) / (S_20 - S_21);
	double delta_11 = (V_21 - V_22) / (S_21 - S_22);

	double theta = (V_21 - V) / (2 * dt);

	double gamma = 2 * (delta_10 - delta_11) / (S_20 - S_22);

	return vector<double>{V, delta, theta, gamma};
}

// using average Binomial Tree to price euro and american options 
vector<double> BT_avg(int N, double S, double K, double r, double T, 
	double q, double sig, const char option, const std::string style)
{
	vector<double> results_N1(4);
	vector<double> results_N2(4);
	if (style == "am")
	{
		results_N1 = BT_am(N, S, K, r, T, q, sig, option);
		results_N2 = BT_am(N + 1, S, K, r, T, q, sig, option);
	}
	else
	{
		results_N1 = BT_euro(N, S, K, r, T, q, sig, option);
		results_N2 = BT_euro(N + 1, S, K, r, T, q, sig, option);
	}
	vector<double>results(4);
	for (int i = 0; i < 4; i++)
	{
		results[i] = (results_N1[i] + results_N2[i]) / 2.0;
	}
	return results;
}

// using BBS to price euro options 
vector<double> BBS_euro(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1.0 / u;
	double P_up = (exp((r - q) * dt) - d) / (u - d);
	double P_down = 1 - P_up;


	vector<double> nodes;
	double d_square = d * d;
	double S_temp = S * pow(u, N - 1);
	double disc = exp(-r * dt);

	for (int i = 0; i < N ; i++)
	{

		auto node_temp = BlackScholes(S_temp, K, r, dt, q, sig, option);
		nodes.push_back(node_temp[0]);
		S_temp *= d_square;
	}

	for (int i = N - 2; i > 1; i--)
	{
		for (int j = 0; j < i + 1; j++)
		{
			nodes[j] = (P_up * nodes[j] + P_down * nodes[j + 1]);
		}
	}

	double S_10 = S * u;
	double S_11 = S * d;

	double S_20 = S * u * u;
	double S_21 = S;
	double S_22 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];

	double V_10 = (P_up * V_20 + P_down * V_21);
	double V_11 = (P_up * V_21 + P_down * V_22);

	double V = exp(-r * T) * (P_up * V_10 + P_down * V_11);

	double delta = exp(-r * T / N * (N - 1)) * (V_10 - V_11) / (S_10 - S_11);
	double delta_10 = (V_20 - V_21) / (S_20 - S_21);
	double delta_11 = (V_21 - V_22) / (S_21 - S_22);

	double theta = (V_21 * exp(-r * T / (double)N * (double)(N - 2)) - V) / (2 * dt);

	double gamma = 2 * exp(-r * T / N * (N - 1)) * (delta_10 - delta_11) / (S_20 - S_22);

	return vector<double>{V, delta, theta, gamma};
}

// using BBS to price am options 
vector<double> BBS_am(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / N;
	double u = exp(sig * sqrt(dt));
	double d = 1.0 / u;
	double P_up = (exp((r - q) * dt) - d) / (u - d);
	double P_down = 1 - P_up;


	vector<double> nodes;
	vector<double> S_list;
	double d_square = d * d;
	double S_temp = S * pow(u, N);
	double disc = exp(-r * dt);

	for (int i = 0; i < N + 1; i++)
	{
		double S_T = S_temp;
		S_list.push_back(S_T);
		S_temp *= d_square;
	}


	for (int i = 0; i < N; i++)
	{
		S_list[i] = S_list[i] * d;
		auto node_temp = BlackScholes(S_list[i], K, r, dt, q, sig, option);
		if (option == 'C')
		{
			nodes.push_back(max(node_temp[0], max(S_list[i]-K, 0)));
		}
		else
		{
			nodes.push_back(max(node_temp[0], max(K- S_list[i], 0)));
		}
	}

	for (int i = N - 2; i > 1; i--)
	{
		for (int j = 0; j < i + 1; j++)
		{
			S_list[j] = S_list[j] * d;
			auto node_temp1 = disc * (P_up * nodes[j] + P_down * nodes[j + 1]);
			double node_temp2 = 0;
			if (option == 'C')
			{
				node_temp2 = max(S_list[j] - K, 0);
			}
			else
			{
				node_temp2 = max(K - S_list[j], 0);
			}

			nodes[j] = max(node_temp1, node_temp2);
		}

	}

	double S_10 = S * u;
	double S_11 = S * d;

	double S_20 = S * u * u;
	double S_21 = S;
	double S_22 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];

	double V_10am = 0;
	double V_11am = 0;
	if (option == 'C')
	{
		V_10am = max(S_10 - K, 0);
		V_11am = max(S_11 - K, 0);
	}
	else
	{
		V_10am = max(K - S_10, 0);
		V_11am = max(K - S_11, 0);
	}
	double V_10 = max(disc * (P_up * V_20 + P_down * V_21), V_10am);
	double V_11 = max(disc * (P_up * V_21 + P_down * V_22), V_10am);

	double V_am = 0;
	if (option == 'C')
	{
		V_am = max(S - K, 0);
	}
	else
	{
		V_am = max(K - S, 0);
	}
	double V = max(disc * (P_up * V_10 + P_down * V_11), V_am);

	double delta = (V_10 - V_11) / (S_10 - S_11);
	double delta_10 = (V_20 - V_21) / (S_20 - S_21);
	double delta_11 = (V_21 - V_22) / (S_21 - S_22);

	double theta = (V_21 - V) / (2 * dt);

	double gamma = 2 * (delta_10 - delta_11) / (S_20 - S_22);

	return vector<double>{V, delta, theta, gamma};
}

// using BBSR to price euro and american options 
vector<double> BBSR(int N, double S, double K, double r, double T, double q,
	double sig, const char option, const std::string style)
{
	vector<double> results_N1(4);
	vector<double> results_N2(4);
	if (style == "am")
	{
		results_N1 = BBS_am(N, S, K, r, T, q, sig, option);
		results_N2 = BBS_am(N / 2, S, K, r, T, q, sig, option);
	}
	else
	{
		results_N1 = BBS_euro(N, S, K, r, T, q, sig, option);
		results_N2 = BBS_euro(N / 2, S, K, r, T, q, sig, option);
	}
	vector<double>results(4);
	for (int i = 0; i < 4; i++)
	{
		results[i] = 2 * results_N1[i] - results_N2[i];
	}
	return results;
}


// using Trinomial Tree to price euro options 
vector<double> TT_euro(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / (double)N;
	double u = exp(sig * sqrt(3*dt));
	double d = 1.0 / u;
	double disc = exp(-r * dt);

	double P_up = (1.0 / 6.0 + (r - q - sig * sig / 2.0) * sqrt(dt / 12.0 / sig / sig)) * disc;
	double P_m = 2.0 / 3.0 * disc;
	double P_down = disc - P_up - P_m;
	

	vector<double> nodes(2 * N + 1);
	double S_up = S;
	double S_down = S;
	if (option == 'C')
	{
		nodes[N] = max(S - K, 0);
	}
	else
	{
		nodes[N] = max(K-S, 0);
	}
	

	for (int i = 1; i < N + 1; i++)
	{
		S_up *= u;
		S_down *= d;
		
		if (option == 'C')
		{
			nodes[N - i] = max(S_up - K, 0);
			nodes[N + i] = max(S_down - K, 0);
		}
		else
		{
			nodes[N - i] = max(K - S_up, 0);
			nodes[N + i] = max(K - S_down, 0);
		}
	}

	for (int i = N - 1; i > 1; i--)
	{
		for (int j = 0; j < 2 * i + 1; j++)
		{
			nodes[j] = (P_up * nodes[j] + P_m * nodes[j + 1] + P_down * nodes[j + 2]);
		}
	}

	
	double S_10 = S * u;
	double S_11 = S;
	double S_12 = S * d;

	double S_20 = S * u * u;
	double S_21 = S * u;
	double S_22 = S;
	double S_23 = S * d;
	double S_24 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];
	double V_23 = nodes[3];
	double V_24 = nodes[4];

	double V_10 = (P_up * V_20 + P_down * V_22 + P_m * V_21);
	double V_11 = (P_up * V_21 + P_down * V_23 + P_m * V_22);
	double V_12 = (P_up * V_22 + P_down * V_24 + P_m * V_23);

	double V =  (P_up * V_10 + P_down * V_12 + P_m * V_11);

	double delta = (V_10 - V_12) / (S_10 - S_12);
	double delta1 = (V_20 - V_22) / (S_20 - S_22);
	double delta2 = (V_22 - V_24) / (S_22 - S_24);

	double theta = (V_11 - V) / dt;

	double gamma = (delta1 - delta2) / (S_20 - S_24) * 2.0;

	return vector<double>{V, delta, theta, gamma};
}

// using Trinomial Black-Scholes Tree to price euro options 
vector<double> TBS_euro(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / (double)N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1.0 / u;
	double disc = exp(-r * dt);

	double P_up = (1.0 / 6.0 + (r - q - sig * sig / 2.0) * sqrt(dt / 12.0 / sig / sig)) * disc;
	double P_m = 2.0 / 3.0 * disc;
	double P_down = disc - P_up - P_m;


	vector<double> nodes(2 * N + 1);
	double S_up = S;
	double S_down = S;
	nodes[N-1] = BlackScholes(S, K, r, dt, q, sig, option)[0];



	for (int i = 1; i < N; i++)
	{
		S_up *= u;
		S_down *= d;
		nodes[N - 1 - i] = BlackScholes(S_up, K, r, dt, q, sig, option)[0];
		nodes[N - 1 + i] = BlackScholes(S_down, K, r, dt, q, sig, option)[0];
	}

	for (int i = N - 2; i > 1; i--)
	{
		for (int j = 0; j < 2 * i + 1; j++)
		{
			nodes[j] = (P_up * nodes[j] + P_m * nodes[j + 1] + P_down * nodes[j + 2]);
		}
	}


	double S_10 = S * u;
	double S_11 = S;
	double S_12 = S * d;

	double S_20 = S * u * u;
	double S_21 = S * u;
	double S_22 = S;
	double S_23 = S * d;
	double S_24 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];
	double V_23 = nodes[3];
	double V_24 = nodes[4];

	double V_10 = (P_up * V_20 + P_down * V_22 + P_m * V_21);
	double V_11 = (P_up * V_21 + P_down * V_23 + P_m * V_22);
	double V_12 = (P_up * V_22 + P_down * V_24 + P_m * V_23);

	double V = (P_up * V_10 + P_down * V_12 + P_m * V_11);

	double delta = (V_10 - V_12) / (S_10 - S_12);
	double delta1 = (V_20 - V_22) / (S_20 - S_22);
	double delta2 = (V_22 - V_24) / (S_22 - S_24);

	double theta = (V_11 - V) / dt;

	double gamma = (delta1 - delta2) / (S_20 - S_24) * 2.0;

	return vector<double>{V, delta, theta, gamma};
}

// using Trinomial Black-Scholes Tree with richard extrapolation to price euro options 
vector<double> TBSR(int N, double S, double K, double r, double T, double q,
	double sig, const char option, const std::string style)
{
	vector<double> results_N1(4);
	vector<double> results_N2(4);
	if (style == "am")
	{
		results_N1 = TBS_am(N, S, K, r, T, q, sig, option);
		results_N2 = TBS_am(N / 2, S, K, r, T, q, sig, option);
	}
	else
	{
		results_N1 = TBS_euro(N, S, K, r, T, q, sig, option);
		results_N2 = TBS_euro(N / 2, S, K, r, T, q, sig, option);
	}
	vector<double>results(4);
	for (int i = 0; i < 4; i++)
	{
		results[i] = 2 * results_N1[i] - results_N2[i];
	}
	return results;
}


// using Trinomial Tree to price am options 
vector<double> TT_am(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / (double)N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1.0 / u;
	double disc = exp(-r * dt);

	double P_up = (1.0 / 6.0 + (r - q - sig * sig / 2.0) * sqrt(dt / 12.0 / sig / sig)) * disc;
	double P_m = 2.0 / 3.0 * disc;
	double P_down = disc - P_up - P_m;


	vector<double> nodes(2 * N + 1);
	vector<double> s_list(2 * N + 1);
	s_list[N] = S;
	double S_up = S;
	double S_down = S;
	if (option == 'C')
	{
		nodes[N] = max(S - K, 0);
	}
	else
	{
		nodes[N] = max(K - S, 0);
	}



	for (int i = 1; i < N + 1; i++)
	{
		S_up *= u;
		S_down *= d;
		s_list[N - i] = S_up;
		s_list[N + i] = S_down;
		if (option == 'C')
		{
			nodes[N - i] = max(S_up - K, 0);
			nodes[N + i] = max(S_down - K, 0);
		}
		else
		{
			nodes[N - i] = max(K - S_up, 0);
			nodes[N + i] = max(K - S_down, 0);
		}
	}

	for (int i = N - 1; i > 1; i--)
	{
		for (int j = 0; j < 2 * i + 1; j++)
		{
			s_list[j] = s_list[j] * d;
			double node_temp = (P_up * nodes[j] + P_m * nodes[j + 1] + P_down * nodes[j + 2]);
			double node_temp2 = 0;
			if (option == 'C')
			{
				node_temp2 = max(s_list[j] - K, 0);
			}
			else
			{
				node_temp2 = max(K - s_list[j], 0);
			}

			nodes[j] = max(node_temp, node_temp2);
		}
	}


	double S_10 = S * u;
	double S_11 = S;
	double S_12 = S * d;

	double S_20 = S * u * u;
	double S_21 = S * u;
	double S_22 = S;
	double S_23 = S * d;
	double S_24 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];
	double V_23 = nodes[3];
	double V_24 = nodes[4];

	double V_10 = (P_up * V_20 + P_down * V_22 + P_m * V_21);
	double V_11 = (P_up * V_21 + P_down * V_23 + P_m * V_22);
	double V_12 = (P_up * V_22 + P_down * V_24 + P_m * V_23);

	double V = (P_up * V_10 + P_down * V_12 + P_m * V_11);

	double delta = (V_10 - V_12) / (S_10 - S_12);
	double delta1 = (V_20 - V_22) / (S_20 - S_22);
	double delta2 = (V_22 - V_24) / (S_22 - S_24);

	double theta = (V_11 - V) / dt;

	double gamma = (delta1 - delta2) / (S_20 - S_24) * 2.0;

	return vector<double>{V, delta, theta, gamma};
}

// using average Trinomial Tree to price euro and american options 
vector<double> TT_avg(int N, double S, double K, double r, double T, double q,
	double sig, const char option, const std::string style)
{
	vector<double> results_N1(4);
	vector<double> results_N2(4);
	if (style == "am")
	{
		results_N1 = TT_am(N, S, K, r, T, q, sig, option);
		results_N2 = TT_am(N + 1, S, K, r, T, q, sig, option);
	}
	else
	{
		results_N1 = TT_euro(N, S, K, r, T, q, sig, option);
		results_N2 = TT_euro(N + 1, S, K, r, T, q, sig, option);
	}
	vector<double>results(4);
	for (int i = 0; i < 4; i++)
	{
		results[i] = (results_N1[i] + results_N2[i]) / 2.0;
	}
	return results;
}


// using Trinomial Black-Scholes Tree to price am options 
vector<double> TBS_am(int N, double S, double K, double r, double T, double q, double sig, const char option)
{
	double dt = (double)T / (double)N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1.0 / u;
	double disc = exp(-r * dt);

	double P_up = (1.0 / 6.0 + (r - q - sig * sig / 2.0) * sqrt(dt / 12.0 / sig / sig)) * disc;
	double P_m = 2.0 / 3.0 * disc;
	double P_down = disc - P_up - P_m;

	vector<double> nodes(2 * N + 1);
	vector<double> s_list(2 * N + 1);
	s_list[N] = S;
	double S_up = S;
	double S_down = S; 

	for (int i = 1; i < N + 1; i++)
	{
		S_up *= u;
		S_down *= d;
		s_list[N - i] = S_up;
		s_list[N + i] = S_down;
	}
	
	
	for (int i = 0; i < 2 * N; i++)
	{
		s_list[i] = s_list[i + 1];
	}

	double node_temp = BlackScholes(S, K, r, dt, q, sig, option)[0];
	double node_temp2 = 0;
	if (option == 'C')
	{
		node_temp2 = max(s_list[N - 1] - K, 0);
	}
	else
	{
		node_temp2 = max(K - s_list[N - 1], 0);
	}
	nodes[N - 1] = max(node_temp, node_temp2);

	double node_temp3 = 0;
	double node_temp4 = 0;
	for (int i = 1; i < N; i++)
	{
		
		node_temp= BlackScholes(s_list[N - 1 - i], K, r, dt, q, sig, option)[0];
		node_temp3= BlackScholes(s_list[N - 1 + i], K, r, dt, q, sig, option)[0];
		if (option == 'C')
		{
			node_temp2 = max(s_list[N - 1 - i] - K, 0);
			node_temp4 = max(s_list[N - 1 + i] - K, 0);
		}
		else
		{
			node_temp2 = max(K - s_list[N - 1 - i], 0);
			node_temp4 = max(K - s_list[N - 1 + i], 0);
		}
		nodes[N - 1 - i] = max(node_temp, node_temp2);
		nodes[N - 1 + i] = max(node_temp3, node_temp4);
	}

	for (int i = N - 2; i > 1; i--)
	{
		for (int j = 0; j < 2 * i + 1; j++)
		{
			s_list[j] = s_list[j + 1];
			node_temp = (P_up * nodes[j] + P_m * nodes[j + 1] + P_down * nodes[j + 2]);
			if (option == 'C')
			{
				node_temp2 = max(s_list[j] - K, 0);
			}
			else
			{
				node_temp2 = max(K - s_list[j], 0);
			}

			nodes[j] = max(node_temp, node_temp2);
		}
	}


	double S_10 = S * u;
	double S_11 = S;
	double S_12 = S * d;

	double S_20 = S * u * u;
	double S_21 = S * u;
	double S_22 = S;
	double S_23 = S * d;
	double S_24 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];
	double V_23 = nodes[3];
	double V_24 = nodes[4];

	double V_10 = (P_up * V_20 + P_down * V_22 + P_m * V_21);
	double V_11 = (P_up * V_21 + P_down * V_23 + P_m * V_22);
	double V_12 = (P_up * V_22 + P_down * V_24 + P_m * V_23);

	double V = (P_up * V_10 + P_down * V_12 + P_m * V_11);

	double delta = (V_10 - V_12) / (S_10 - S_12);
	double delta1 = (V_20 - V_22) / (S_20 - S_22);
	double delta2 = (V_22 - V_24) / (S_22 - S_24);

	double theta = (V_11 - V) / dt;

	double gamma = (delta1 - delta2) / (S_20 - S_24) * 2.0;

	return vector<double>{V, delta, theta, gamma};
}

// using trinomial tree to price DaO calls 
vector<double> TT_euro_DAOC(int N, double S, double K, double r, double T, double q, double sig, double B, const char option)
{
	double dt = (double)T / (double)N;
	double u = exp(sig * sqrt(3 * dt));
	double d = 1.0 / u;
	double disc = exp(-r * dt);

	double P_up = (1.0 / 6.0 + (r - q - sig * sig / 2.0) * sqrt(dt / 12.0 / sig / sig)) * disc;
	double P_m = 2.0 / 3.0 * disc;
	double P_down = disc - P_up - P_m;


	vector<double> nodes(2 * N + 1);
	vector<double> s_list(2 * N + 1);
	s_list[N] = S;
	double S_up = S;
	double S_down = S;
	nodes[N] = max(S - K, 0);



	for (int i = 1; i < N + 1; i++)
	{
		S_up *= u;
		S_down *= d;
		s_list[N - i] = S_up;
		s_list[N + i] = S_down;
		nodes[N - i] = max(S_up - K, 0);
		nodes[N + i] = max(S_down - K, 0);
	}

	for (int i = N - 1; i > 1; i--)
	{
		for (int j = 0; j < 2 * i + 1; j++)
		{
			s_list[j] = s_list[j] * d;
			double node_temp = (P_up * nodes[j] + P_m * nodes[j + 1] + P_down * nodes[j + 2]);
			if (s_list[j] > B)
			{
				nodes[j] = node_temp;
			}
			else
			{
				nodes[j] = 0;
			}

			
		}
	}

	double S_10 = S * u;
	double S_11 = S;
	double S_12 = S * d;

	double S_20 = S * u * u;
	double S_21 = S * u;
	double S_22 = S;
	double S_23 = S * d;
	double S_24 = S * d * d;

	double V_20 = nodes[0];
	double V_21 = nodes[1];
	double V_22 = nodes[2];
	double V_23 = nodes[3];
	double V_24 = nodes[4];

	double V_10 = (P_up * V_20 + P_down * V_22 + P_m * V_21);
	double V_11 = (P_up * V_21 + P_down * V_23 + P_m * V_22);
	double V_12 = (P_up * V_22 + P_down * V_24 + P_m * V_23);

	double V = (P_up * V_10 + P_down * V_12 + P_m * V_11);

	double delta = (V_10 - V_12) / (S_10 - S_12);
	double delta1 = (V_20 - V_22) / (S_20 - S_22);
	double delta2 = (V_22 - V_24) / (S_22 - S_24);

	double theta = (V_11 - V) / dt;

	double gamma = (delta1 - delta2) / (S_20 - S_24) * 2.0;


	return vector<double>{V, delta, theta, gamma};
}


