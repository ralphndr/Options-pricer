#include "../include/AmericanSolver.hpp"
#include <vector>
#include <algorithm>

AmericanSolver::AmericanSolver(const AmericanOption& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

double AmericanSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dt = grid.getdt();
    double dS = grid.getdS();

    // Initialisation à maturité
    for(int i = 0; i <= M; i++)
        grid.set(i, N, option.payoff(i * dS));

    std::vector<double> a(M-1), b(M-1), c(M-1), d(M-1);

    // Schéma implicite avec exercice anticipé
    for(int j = N-1; j >= 0; j--) {
        for(int i = 1; i < M; i++) {
            double S = i * dS;
            double alpha = 0.25 * dt * (data.sigma*data.sigma*S*S/dS/dS - data.r*S/dS);
            double beta  = -0.5 * dt * (data.sigma*data.sigma*S*S/dS/dS + data.r);
            double gamma = 0.25 * dt * (data.sigma*data.sigma*S*S/dS/dS + data.r*S/dS);

            a[i-1] = -alpha;
            b[i-1] = 1 - beta;
            c[i-1] = -gamma;

            d[i-1] = alpha*grid.get(i-1,j+1) + (1+beta)*grid.get(i,j+1) + gamma*grid.get(i+1,j+1);
        }

        // Résolution tridiagonale (Thomas)
        std::vector<double> c_star(M-1), d_star(M-1);
        c_star[0] = c[0]/b[0];
        d_star[0] = d[0]/b[0];
        for(int i = 1; i < M-2; i++) {
            double m = b[i] - a[i]*c_star[i-1];
            c_star[i] = c[i]/m;
            d_star[i] = (d[i] - a[i]*d_star[i-1])/m;
        }
        d[M-2] = d_star[M-2];
        for(int i = M-3; i >= 0; i--) {
            d[i] = d_star[i] - c_star[i]*d[i+1];
        }

        // Mise à jour avec exercice anticipé
        for(int i = 1; i < M; i++)
            grid.set(i, j, std::max(d[i-1], option.payoff(i*dS)));

        grid.set(0, j, option.payoff(0.0));
        grid.set(M, j, option.payoff(M*dS));
    }

    int iStrike = static_cast<int>(option.getStrike()/dS);
    return grid.get(iStrike, 0);
}
