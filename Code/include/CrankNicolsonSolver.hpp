#ifndef CRANK_NICOLSON_SOLVER_HPP
#define CRANK_NICOLSON_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

class CrankNicolsonSolver : public PDESolver {
public:
    CrankNicolsonSolver(const Option& opt, const MarketData& market,
                        int M, int N, double Smax);

    double price() override; // doit être défini dans le .cpp

private:
    void setupTridiagonalCN(std::vector<double>& a,
                             std::vector<double>& b,
                             std::vector<double>& c,
                             double dt);
    void solveTridiagonal(const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& c,
                          std::vector<double>& d);
};

#endif
