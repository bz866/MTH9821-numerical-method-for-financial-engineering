#pragma once
#include"FiniteDifferencesPricer.h"

void problem2() {
	double S0, K, T, sigma, q, r;
	S0 = 42.;
	K = 40.;
	T = 0.75;
	sigma = 0.32;
	q = 0.02;
	r = 0.04;

	Option op(American, put, S0, K, T, sigma, q, r);
	FiniteDifferencePutOptionPricer pricer_1(op, 0.45, 4, Forward);
	FiniteDifferencePutOptionPricer pricer_2(op, 0.45, 4, CN);
	
	ofstream out;
	out.open("2_data.csv");
	out << "AMERICAN" << endl; 

	for (auto & pricer : { pricer_1,pricer_2 }) {
		if (pricer.fd_type == 0)out << "Forward" << endl;
		else out << "CN" << endl;
		for (int m = 0;m <= pricer.M;++m) {
			for (int n = 0;n <= pricer.N;++n) {
				out << setprecision(10) << pricer.Domin[n][m] << ",";
			}
			out << endl;
		}
		out << endl << endl;
	}

	// the first part: output domin end
	///////////////////////////////////////////
	

	for (double alpha : {0.45, 5.0}) {
		for (auto fd_t : { Forward, CN}) {
			out << alpha << ",";
			if (fd_t == 0)out << "Forward" << endl;
			else out << "CN" << endl;
			for (int M_ = 4;M_ <= 256;M_ *= 4) {
				FiniteDifferencePutOptionPricer p(op, alpha, M_, fd_t);
				p.computeRMS();
				p.computeGreeks();
				out << setprecision(10) << p.error_pointwise << "," << "," <<
					p.error_pointwise2 << "," << "," <<
					p.RMS << "," << "," <<
					p.delta_approx << "," << p.gamma_approx << "," << p.theta_approx << endl;
			}
			out << endl << endl << endl;
		}
	}
	
	
	out.close();
	
	return;
}
