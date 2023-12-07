#ifndef TREEHPP
#define TREE_HPP

#include<vector>
#include<string>
#include <iostream>

using namespace std;

// return max of two numbers 
double max(double d1, double d2);

// using Binomial Tree to price euro options 
vector<double> BT_euro(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using Binomial Tree to price american options 
vector<double> BT_am(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using average Binomial Tree to price euro and american options 
vector<double> BT_avg(int N, double S, double K, double r, double T, double q, 
	double sig, const char option, const std::string style);

// using BBS to price euro options 
vector<double> BBS_euro(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using BBS to price am options 
vector<double> BBS_am(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using BBSR to price euro and american options 
vector<double> BBSR(int N, double S, double K, double r, double T, double q,
	double sig, const char option, const std::string style);

// using Trinomial Tree to price euro options 
vector<double> TT_euro(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using Trinomial Black-Scholes Tree to price euro options 
vector<double> TBS_euro(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using Trinomial Black-Scholes Tree with richard extrapolation to price euro and am options 
vector<double> TBSR(int N, double S, double K, double r, double T, double q, 
	double sig, const char option, const std::string style);

// using Trinomial Tree to price am options 
vector<double> TT_am(int N, double S, double K, double r, double T, double q, double sig, const char option);

// using Trinomial Black-Scholes Tree to price am options 
vector<double> TBS_am(int N, double S, double K, double r, double T, double q, double sig, const char option);


// using average Trinomial Tree to price euro and american options 
vector<double> TT_avg(int N, double S, double K, double r, double T, double q,
	double sig, const char option, const std::string style);

// using trinomial tree to price DaO calls 
vector<double> TT_euro_DAOC(int N, double S, double K, double r, double T, double q, double sig, double B, const char option);



#endif // !TREE_HPP
