#include "../include/PDESolver.hpp"

PDESolver::PDESolver(const Option& opt,
                     const MarketData& market,
                     int M, int N, double Smax)
    : option(opt), data(market), grid(M, N, Smax, opt.getMaturity()) {}
