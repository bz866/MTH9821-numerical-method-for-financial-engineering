#ifndef UTILS_HPP
#define UTILS_HPP


#include<iostream>
#include<vector>

inline void printVector(const std::vector<double> &vec)
{
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        std::cout << vec[i] << "\t";
    }
    std::cout << std::endl;
}

#endif
