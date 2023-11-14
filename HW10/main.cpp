#include "EuropeanOption.hpp"
#include "FiniteDifferencePricer.hpp"
#include <iostream>
#include <iomanip>


int main()
{
    // ================  Question 1 ===========================
    // Pricing European Put Options Using Finite Differences on a Fixed Computational Domain
    // ========================================================

    std::cout << std::setprecision(8);

    double S0 = 36.;
    double K = 40;
    double T = 0.75;
    double sigma = 0.32;
    double r = 0.04;
    double q = 0.015;

    double alpha = 0.45;
    int M = 4;

    EuropeanOption europeanOption(0, S0, K, T, sigma, r, q);
    FiniteDifferencePricer pricerEuroPut(S0, K, T, sigma, r, q);
    // EuropeanOption europeanOption(0., 37., 40, .75, .28, .03, .015);
    // FiniteDifferencePricer pricerEuroPut(37., 40, .75, .28, .03, .015);
    pricerEuroPut.setHyperparameters(alpha, M);
    
    std::cout << "European Put" << std::endl;
    std::cout << "BS Price: " << europeanOption.Put() << std::endl;

    for (std::size_t i=4; i < 257; i <<=2)
    {
        printVector(pricerEuroPut.priceEuropeanPut(i));
    }


    // ================  Question 2 ===========================
    // Pricing American Put Options Using Finite Differences on a Fixed Computational Domain
    // ========================================================
    S0 = 42;
    K = 40;
    T = 0.75;
    sigma = 0.32;
    r = 0.04;
    q = 0.02;

    EuropeanOption option(0., S0, K, T, sigma, r, q);
    FiniteDifferencePricer pricerAmeriPut(S0, K, T, sigma, r, q);
    pricerAmeriPut.setHyperparameters(alpha, M);
    
    std::cout << "\n\nAmerican Put" << std::endl;
    std::cout << "BS Price: " << 3.3045362802172642 << std::endl;
    
    for (std::size_t i = 4; i < 257; i <<= 2) {
        printVector(pricerAmeriPut.priceAmericanPut(i));
    }


    // Early Exercise
    FiniteDifferencePricer pricerAmericanPutEarlyExercise(S0, K, T, sigma, r, q);
    pricerAmericanPutEarlyExercise.setHyperparameters(alpha, M);
    
    std::cout << "\n\nAmerican Put Early Exercise" << std::endl;
    
    std::vector<double> t;
    std::vector<double> Sopt;
    std::tie(t, Sopt) = pricerAmericanPutEarlyExercise.priceAmericanPutEarlyExercise(16);
    printVector(t);
    printVector(Sopt);


    // // ================  Question 3 ===========================
    // // Pricing American Put Options Using Finite Differences on a Fixed Computational Domain
    // // ========================================================
    // S0 = 36;
    // K = 40;
    // T = 0.75;
    // sigma = 0.32;
    // r = 0.045;
    // q = 0.02;

    // EuropeanOption optionQ3(0., S0, K, T, sigma, r, q);
    // FiniteDifferencePricer priceEuroPut_BE(S0, K, T, sigma, r, q);
    // priceEuroPut_BE.setHyperparameters(alpha, M);

    // std::cout << "\n\nEuropean Put Backward Euler" << std::endl;
    // std::cout << "BS Price: " << optionQ3.Put() << std::endl;
    // for (std::size_t i=4; i < 257; i <<=2)
    // {
    //     printVector(priceEuroPut_BE.priceEuropeanPut_BE(i));
    // }

    // // ================  Question 4 ===========================
    // // Pricing American Put Options Using Finite Differences on a Fixed Computational Domain
    // // ========================================================
    // S0 = 36;
    // K = 40;
    // T = 0.75;
    // sigma = 0.32;
    // r = 0.045;
    // q = 0.02;

    // EuropeanOption optionQ4(0., S0, K, T, sigma, r, q);
    // FiniteDifferencePricer pricerAmeriPut_BE(S0, K, T, sigma, r, q);
    // pricerAmeriPut_BE.setHyperparameters(alpha, M);

    // std::cout << "\n\nAmerican Put Backward Euler" << std::endl;
    // std::cout << "BS Price: " << optionQ4.Put() << std::endl;
    // for (std::size_t i=4; i < 257; i <<=2)
    // {
    //     printVector(pricerAmeriPut_BE.priceEuropeanPut_BE(i));
    // }

    // HW10 
    // ================================= Question 1 =================================
    S0 = 37;
    K = 40;
    T = 0.75;
    sigma = 0.28;
    r = 0.04;
    q = 0.02;
    EuropeanOption option_HW10_Q1(0., S0, K, T, sigma, r, q);
    FiniteDifferencePricer priceAmeriPut_
    
    return 0;
}
