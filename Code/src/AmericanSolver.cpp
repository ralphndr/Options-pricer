#include "../include/AmericanSolver.hpp"
#include <vector>
#include <algorithm>

// Constructeur du solveur pour option américaine.
// Initialise la grille et les paramètres via le constructeur de PDESolver.
AmericanSolver::AmericanSolver(const AmericanOption& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

// Calcule le prix de l'option américaine via un schéma implicite
// avec gestion de l'exercice anticipé (condition d'obstacle).
double AmericanSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dt = grid.getdt();
    double dS = grid.getdS();
    double r = data.r;
    double T = N * dt;

    // Initialisation à maturité : payoff de l'option
    for(int i = 0; i <= M; i++)
        grid.set(i, N, option.payoff(i * dS));

    std::vector<double> a(M-1), b(M-1), c(M-1), d(M-1);

    // Détermination du type de l'option pour les conditions aux bords
    const AmericanOption& am_opt = static_cast<const AmericanOption&>(option);
    bool is_call = (am_opt.getType() == AmericanOptionType::Call);

    // Boucle arrière en temps (schéma implicite avec exercice anticipé)
    for(int j = N-1; j >= 0; j--) {
        double tau = j * dt;
        double time_to_expiry = T - tau;
        
        // Construction des coefficients de la matrice tridiagonale et du vecteur d
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

        // Résolution tridiagonale par l'algorithme de Thomas
        std::vector<double> c_star(M-1), d_star(M-1);
        c_star[0] = c[0]/b[0];
        d_star[0] = d[0]/b[0];
        for(int i = 1; i < M-1; i++) {
            double m = b[i] - a[i]*c_star[i-1];
            c_star[i] = c[i]/m;
            d_star[i] = (d[i] - a[i]*d_star[i-1])/m;
        }
        d_star[M-2] = (d[M-2] - a[M-2]*d_star[M-3])/(b[M-2] - a[M-2]*c_star[M-3]);
        
        d[M-2] = d_star[M-2];
        for(int i = M-3; i >= 0; i--) {
            d[i] = d_star[i] - c_star[i]*d[i+1];
        }

        // Mise à jour de la grille avec condition d'exercice anticipé
        for(int i = 1; i < M; i++)
            grid.set(i, j, std::max(d[i-1], option.payoff(i*dS)));

        // Conditions aux bords (mêmes formules asymptotiques que pour l'option européenne)
        double S_max = M * dS;
        if(is_call) {
            // Call : V(0,t) = 0, V(Smax,t) = Smax - K*exp(-r*(T-t))
            grid.set(0, j, 0.0);
            grid.set(M, j, S_max - option.getStrike() * std::exp(-r * time_to_expiry));
        } else {
            // Put : V(0,t) = K*exp(-r*(T-t)), V(Smax,t) = 0
            grid.set(0, j, option.getStrike() * std::exp(-r * time_to_expiry));
            grid.set(M, j, 0.0);
        }
    }

    // Interpolation grossière pour récupérer la valeur au strike
    int iStrike = static_cast<int>(option.getStrike()/dS);
    return grid.get(iStrike, 0);
}
