#include "NumericalMethod.hpp"

#include <cmath>

using namespace std;

// Linear Congruential Generator
vector<long double> NumericalMethod::LinearCongruentialRandom(unsigned long int size) {
    vector<long double> rand(size);
    
    unsigned long int a = 39373;
    unsigned long int c = 0;
    unsigned long int k = pow(2, 31) - 1;
    unsigned long int x0 = 1;
    unsigned long int x = x0;
    for (unsigned long int i = 0; i < size; ++i) {
        x = (a * x + c) % k;
        rand[i] = (long double) x / (long double) k;
    }
    return rand;
}

// Inverse Transform Method
vector<long double> NumericalMethod::InverseTransformRandom(unsigned long int size) {
    vector<long double> rand(size);
    vector<long double> u = LinearCongruentialRandom(size);
    long double x, y, r;
    
    double a0 = 2.50662823884;
    double a1 = -18.61500062529;
    double a2 = 41.39119773534;
    double a3 = -25.44106049637;
    double b0 = -8.47351093090;
    double b1 = 23.08336743743;
    double b2 = -21.06224101826;
    double b3 = 3.13082909833;
    double c0 = 0.3374754822726147;
    double c1 = 0.9761690190917186;
    double c2 = 0.16079797149118209;
    double c3 = 0.0276438810333863;
    double c4 = 0.0038405729373609;
    double c5 = 0.0003951896511919;
    double c6 = 0.0000321767881768;
    double c7 = 0.0000002888167364;
    double c8 = 0.0000003960315187;

    for (unsigned long int i = 0; i < size; ++i) {
        y = u[i] - 0.5;
        if (abs(y) < 0.42) {
            r = y * y;
            x = y * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.0);
        } else {
            r = u[i];
            if (y > 0) { r = 1 - u[i]; }
            r = log(-log(r));
            x = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))));
            if (y < 0) { x = -x; }
        }
        rand[i] = x;
    }
    
    return rand;
}

// Acceptance-Rejection Method
vector<long double> NumericalMethod::AcceptanceRejectionRandom(unsigned long int size) {
    vector<long double> rand(size);
    vector<long double> u = LinearCongruentialRandom(size*3);
    long double X;
    unsigned long int n = 0;
    for (unsigned long int i = 0; i < u.size()-2; i += 3) {
        X = -log(u[i]);
        if (u[i+1] > exp(-0.5 * (X-1) * (X-1))) {
            continue;
        } else {
            if (u[i+2] <= 0.5) { X = -X; }
        }
        
        rand[n] = X;
        n++;
    }
    rand.resize(n);
    
    return rand;
}

// Box-Muller Method
vector<long double> NumericalMethod::BoxMullerRandom(unsigned long int size) {
    vector<long double> rand(2*size);
    vector<long double> u = LinearCongruentialRandom(2*size);
    long double u1, u2, X, Y;
    unsigned long int n = 0;
    for (unsigned long int i=0; i < u.size()-1; i += 2) {
        u1 = 2 * u[i] - 1;
        u2 = 2 * u[i+1] - 1;
        X = u1 * u1 + u2 * u2;
        if (X > 1) {
            continue;
        }
        
        Y = sqrt(-2 * log(X) / X);        
        rand[n] = u1 * Y;
        rand[n+1] = u2 * Y;
        n += 2;
    }
    rand.resize(n);

    return rand;
}
