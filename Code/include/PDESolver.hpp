#ifndef PDE_SOLVER_HPP
#define PDE_SOLVER_HPP

#include "Option.hpp"
#include "MarketData.hpp"
#include "Grid.hpp"

class PDESolver {
public:
    PDESolver(const Option& opt, const MarketData& market,
              int M, int N, double Smax);

    virtual ~PDESolver() = default;

    virtual double price() = 0;

protected:
    const Option& option;
    const MarketData& data;
    Grid grid;
};

#endif
