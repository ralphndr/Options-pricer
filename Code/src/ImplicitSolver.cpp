#include "../include/ImplicitSolver.hpp"
#include "../include/EuropeanOption.hpp"
#include <vector>
#include <iostream>

ImplicitSolver::ImplicitSolver(const Option& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

void ImplicitSolver::setupTridiagonal(std::vector<double>& a,
                                      std::vector<double>& b,
                                      std::vector<double>& c,
                                      double dt) {
    int M = grid.getSizeS();
    double dS = grid.getdS();
    double r = data.r;
    double sigma = data.sigma;

    a.resize(M-1);
    b.resize(M-1);
    c.resize(M-1);

    for(int i=1; i<M; ++i){
        double S = i * dS;
        a[i-1] = -0.5*dt*(sigma*sigma*S*S/dS/dS - r*S/(2*dS));
        b[i-1] = 1 + dt*(sigma*sigma*S*S/dS/dS + r);
        c[i-1] = -0.5*dt*(sigma*sigma*S*S/dS/dS + r*S/(2*dS));
    }
}

// Méthode Thomas pour résoudre tridiagonale
void ImplicitSolver::solveTridiagonal(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      const std::vector<double>& c,
                                      std::vector<double>& d){
    int n = b.size();
    std::vector<double> c_star(n);
    std::vector<double> d_star(n);

    c_star[0] = c[0]/b[0];
    d_star[0] = d[0]/b[0];

    for(int i=1; i<n; ++i){
        double m = b[i] - a[i]*c_star[i-1];
        c_star[i] = c[i]/m;
        d_star[i] = (d[i] - a[i]*d_star[i-1])/m;
    }

    d[n-1] = d_star[n-1];
    for(int i=n-2;i>=0;i--){
        d[i] = d_star[i] - c_star[i]*d[i+1];
    }
}

double ImplicitSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dt = grid.getdt();
    double dS = grid.getdS();

    for(int i=0;i<=M;++i)
        grid.set(i,N,option.payoff(i*dS));

    std::vector<double> a,b,c,d(M-1);

    for(int j=N-1;j>=0;j--){
        setupTridiagonal(a,b,c,dt);

        for(int i=1;i<M;i++) d[i-1] = grid.get(i,j+1);

        solveTridiagonal(a,b,c,d);

        for(int i=1;i<M;i++) grid.set(i,j,d[i-1]);

        grid.set(0,j,option.payoff(0.0));
        grid.set(M,j,option.payoff(M*dS));
    }

    int iStrike = static_cast<int>(option.getStrike()/dS);
    return grid.get(iStrike,0);
}
