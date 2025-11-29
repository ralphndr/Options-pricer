#ifndef IMPLICIT_SOLVER_HPP
#define IMPLICIT_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

class ImplicitSolver : public PDESolver {
public:
    ImplicitSolver(const Option& opt, const MarketData& market,
                   int M, int N, double Smax);

    double price() override;

private:
    void setupTridiagonal(std::vector<double>& a,
                          std::vector<double>& b,
                          std::vector<double>& c,
                          double dt);
    void solveTridiagonal(const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& c,
                          std::vector<double>& d);
};

#endif
