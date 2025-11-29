#include "../include/ExplicitSolver.hpp"
#include "../include/EuropeanOption.hpp"
#include <algorithm>
#include <iostream>

ExplicitSolver::ExplicitSolver(const Option& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

void ExplicitSolver::initializeBoundary() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();

    // Valeurs à maturité (payoff)
    for(int i = 0; i <= M; ++i) {
        double S = i * dS;
        grid.set(i, N, option.payoff(S));
    }

    // Conditions aux bornes
    grid.set(0, N, option.payoff(0.0));
    grid.set(M, N, option.payoff(M * dS));
}

double ExplicitSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();

    double r = data.r;
    double sigma = data.sigma;

    initializeBoundary();

    // Boucle temps arrière
    for(int j = N - 1; j >= 0; --j) {
        for(int i = 1; i < M; ++i) {
            double S = i * dS;
            double delta = (grid.get(i+1,j+1) - grid.get(i-1,j+1)) / (2*dS);
            double gamma = (grid.get(i+1,j+1) - 2*grid.get(i,j+1) + grid.get(i-1,j+1)) / (dS*dS);

            double V = grid.get(i,j+1) + dt * (0.5 * sigma*sigma*S*S*gamma + r*S*delta - r*grid.get(i,j+1));
            grid.set(i,j, V);
        }
        // Bornes
        grid.set(0,j, option.payoff(0.0));
        grid.set(M,j, option.payoff(M*dS));
    }

    // Interpolation pour S = strike
    int iStrike = static_cast<int>(option.getStrike() / dS);
    return grid.get(iStrike, 0);
}
