#include <iostream>
#include <iomanip>
#include <math.h>
#include <tuple>
#include <string>
#include <Dense>
#include "OptionPricing.hpp"
#include "BinomialPricing.hpp"
using namespace std;
using namespace Eigen;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;


void PrintVector(vector<double> Input) {
	for_each(Input.begin(), Input.end(), [](double i) {cout << setprecision(6) << i << "\t"; });
	cout << endl;
}

void BTTest() {
	double S = 50;
	double K = 55.55;
	double sigma = 0.25;
	double r = 0.04;
	double q = 0.015;
	vector<double> times = { 2.0 / 12, 4.0 / 12, 6.0 / 12 };

	vector<double> OLD, NEW;
	int N;
	double tol;


	
	//Part 1
	cout << "Part 1: " << endl;
	Binomial SOpt(S, K, 3.0 / 12, r, q, sigma);
	Binomial LOpt(S, K, 7.0 / 12, r, q, sigma);
	
	//European Prop 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.European_PropDiv("PUT", N);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.European_PropDiv("PUT", N);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	//American Prop 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.American_PropDiv("PUT", N, 2.0/12);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.American_PropDiv("PUT", N, 2.0/12);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	
	
	//European Prop 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.European_PropDiv_Mult("PUT", N, 3);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW =LOpt.European_PropDiv_Mult("PUT", N, 3);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "European: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	//American  Prop 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.American_PropDiv_Mult("PUT", N, times);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_PropDiv_Mult("PUT", N, times);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	
	
	
	//Part2
	cout << "Part 2: " << endl;
	double v_div = 0.95;
	
	//European Fix 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.European_FixDiv("PUT", N, 2., v_div);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.European_FixDiv("PUT", N, 2.0/12, v_div);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	//American Fix 3M
	tol = 10.1;
	N = 6;
	OLD = SOpt.American_FixDiv("PUT", N, 2.0/12, v_div);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = SOpt.American_FixDiv("PUT", N, 2.0/12, v_div);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "American: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	
	//European Fix 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.European_FixDiv_Mult("PUT", N, times, 0.5);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.European_FixDiv_Mult("PUT", N, times, 0.5);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "European: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	//American Fix 7M
	tol = 10.1;
	N = 7;
	OLD = LOpt.American_FixDiv_Mult("PUT", N, times, 0.5);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_FixDiv_Mult("PUT", N, times, 0.5);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;

	
	//part 3
	vector<double> values = { 0.95, 0.02, 0.45 };
	//European Complex 7M
	cout << "Part 3:" << endl;
	N = 7;
	tol = 10.1;
	OLD = LOpt.European_Complex("PUT", N, times, values);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.European_Complex("PUT", N, times, values);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
	}
	cout << "European: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	//American Complex 7M
	N = 7;
	tol = 10.1;
	
	OLD = LOpt.American_Complex("PUT", N, times, values);
	while (tol >= pow(10, -4)) {
		N = 2 * N;
		NEW = LOpt.American_Complex("PUT", N, times, values);
		tol = abs(NEW[0] - OLD[0]);
		OLD = NEW;
		cout << tol << endl;
	}
	cout << "American: " << setprecision(6) << OLD[0] << "\t" << N << "\t" << setprecision(6) << OLD[1]
		<< "\t" << setprecision(6) << OLD[2] << "\t" << setprecision(6) << OLD[3] << endl;
	
}

int main() {
	cout << "Binomial Tree: " << endl;
	BTTest();
	
	return 0;

}