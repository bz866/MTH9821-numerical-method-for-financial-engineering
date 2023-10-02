#ifndef Calculator_hpp
#define Calculator_hpp

#include <iostream>
#include <iomanip>
#include "TrinomialTree.hpp"
#include "BinomialTreeVolatilityUnknown.hpp"
#include "BinomialTree.hpp"
#include "EuropeanOption.hpp"


class Calculator
{
public:
    Calculator(){};
    ~Calculator(){};
    
    void PriceAndGreeks(const EuropeanOption& option, bool american);
    void PriceAndGreeks_compensated(const EuropeanOption& option);
};

#endif // Calculator.hpp