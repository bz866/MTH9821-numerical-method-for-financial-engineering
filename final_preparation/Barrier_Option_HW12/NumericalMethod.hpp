#ifndef NUMERICAL_METHOD_HPP
#define NUMERICAL_METHOD_HPP

#include <iostream>
#include <vector>

class NumericalMethod {
  public:
    NumericalMethod(){};
    virtual ~NumericalMethod(){};

    std::vector<long double> LinearCongruentialRandom(unsigned long int size);
    std::vector<long double> InverseTransformRandom(unsigned long int size);
    std::vector<long double> AcceptanceRejectionRandom(unsigned long int size);
    std::vector<long double> BoxMullerRandom(unsigned long int size);
};

#endif /* NUMERICAL_METHOD_HPP */
