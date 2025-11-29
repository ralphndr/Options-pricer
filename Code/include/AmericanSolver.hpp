#ifndef AMERICANSOLVER_HPP
#define AMERICANSOLVER_HPP

#include "AmericanOption.hpp"
#include "PDESolver.hpp"

class AmericanSolver : public PDESolver {
public:
    AmericanSolver(const AmericanOption& opt, const MarketData& market,
                   int M, int N, double Smax);
    double price() override;
};

#endif 
