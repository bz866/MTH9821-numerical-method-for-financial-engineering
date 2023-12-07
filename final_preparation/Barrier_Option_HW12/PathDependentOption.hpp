#ifndef PATH_DEPENDENT_OPTION_HPP
#define PATH_DEPENDENT_OPTION_HPP

#include "OptionBase.hpp"

#include <iostream>

class PathDependentOption {
  public:
    PathDependentOption(BarrierOption opt) { option_base = opt; };
    virtual ~PathDependentOption() {};
    
    double ClosedFormCdao();
    double MonteCarloCdao(unsigned long int n, unsigned long int m);

  private:
    BarrierOption option_base;
};

#endif /* PATH_DEPENDENT_OPTION_HPP */
