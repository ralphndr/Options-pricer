#ifndef EXPLICIT_SOLVER_HPP
#define EXPLICIT_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

class ExplicitSolver : public PDESolver {
public:
    ExplicitSolver(const Option& opt, const MarketData& market,
                   int M, int N, double Smax);

    double price() override;

private:
    void initializeBoundary();
};

#endif
