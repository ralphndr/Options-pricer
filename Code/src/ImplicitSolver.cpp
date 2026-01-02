#include "../include/ImplicitSolver.hpp"
#include "../include/VanillaOption.hpp"
#include "../include/ExoticOption.hpp"
#include <vector>
#include <iostream>

// Constructeur du solveur implicite.
// Initialise l'option, les données de marché et la grille via le constructeur PDESolver.
ImplicitSolver::ImplicitSolver(const Option& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

// Construction des coefficients tridiagonaux pour le schéma implicite
// a, b, c : coefficients de la matrice tridiagonale
// dt : pas de temps
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
        a[i-1] = -0.5*dt*(sigma*sigma*S*S/dS/dS - r*S/dS);
        b[i-1] = 1 + dt*(sigma*sigma*S*S/dS/dS + r);
        c[i-1] = -0.5*dt*(sigma*sigma*S*S/dS/dS + r*S/dS);
    }
}

// Résolution tridiagonale (méthode de Thomas)
// a, b, c : coefficients de la matrice
// d : vecteur de droite, remplacé par la solution
void ImplicitSolver::solveTridiagonal(const std::vector<double>& a,
                                      const std::vector<double>& b,
                                      const std::vector<double>& c,
                                      std::vector<double>& d){
    int n = b.size();
    std::vector<double> c_star(n);
    std::vector<double> d_star(n);

    // Forward sweep
    c_star[0] = c[0]/b[0];
    d_star[0] = d[0]/b[0];
    for(int i=1; i<n; ++i){
        double m = b[i] - a[i]*c_star[i-1];
        c_star[i] = c[i]/m;
        d_star[i] = (d[i] - a[i]*d_star[i-1])/m;
    }

    // Backward substitution
    d[n-1] = d_star[n-1];
    for(int i=n-2;i>=0;i--){
        d[i] = d_star[i] - c_star[i]*d[i+1];
    }
}

// Calcul du prix de l'option par schéma implicite
double ImplicitSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dt = grid.getdt();
    double dS = grid.getdS();
    double r = data.r;
    double T = N * dt;

    // Initialisation de la grille à maturité avec le payoff
    for(int i=0;i<=M;++i)
        grid.set(i,N,option.payoff(i*dS));

    // Fonction lambda pour appliquer les barrières si option barrière
    auto apply_barrier_row = [&](int j) {
        if (const BarrierOption* bo = dynamic_cast<const BarrierOption*>(&option)) {
            bool knockIn = bo->isKnockIn();
            bool call = (bo->getOptionType() == OptionType::Call);
            double B = bo->getBarrierLevel();
            for (int i = 0; i <= M; ++i) {
                double S = i * dS;
                bool active = knockIn
                    ? (call ? (S >= B) : (S <= B))
                    : (call ? (S <  B) : (S >  B));
                if (!active) grid.set(i, j, 0.0);
            }
        }
    };

    // Application de la barrière à maturité
    apply_barrier_row(N);

    std::vector<double> a,b,c,d(M-1);

    // Boucle arrière en temps
    for(int j=N-1;j>=0;j--){
        // Mise en place des coefficients tridiagonaux
        setupTridiagonal(a,b,c,dt);
        double tau = j * dt;
        double time_to_expiry = T - tau;

        // Détermination du type d'option pour les conditions aux bornes
        bool is_call = true;
        if (const EuropeanOption* euro_opt = dynamic_cast<const EuropeanOption*>(&option)) {
            is_call = euro_opt->isCall();
        } else if (const BarrierOption* bo = dynamic_cast<const BarrierOption*>(&option)) {
            is_call = (bo->getOptionType() == OptionType::Call);
        }
        
        // Conditions aux bornes
        double lower_bc, upper_bc;
        double S_max = M * dS;
        if (const DigitalOption* dig = dynamic_cast<const DigitalOption*>(&option)) {
            bool call = (dig->getOptionType() == OptionType::Call);
            if (call) {
                lower_bc = 0.0;
                upper_bc = dig->getPayout() * std::exp(-r * time_to_expiry);
            } else {
                lower_bc = dig->getPayout() * std::exp(-r * time_to_expiry);
                upper_bc = 0.0;
            }
        } else {
            if(is_call) {
                lower_bc = 0.0;
                upper_bc = S_max - option.getStrike() * std::exp(-r * time_to_expiry);
            } else {
                lower_bc = option.getStrike() * std::exp(-r * time_to_expiry);
                upper_bc = 0.0;
            }
        }

        // Construction du vecteur RHS pour les points intérieurs
        for(int i=1;i<M;i++) {
            d[i-1] = grid.get(i,j+1);
            if(i==1) d[i-1] -= a[i-1] * lower_bc;
            if(i==M-1) d[i-1] -= c[i-1] * upper_bc;
        }

        // Résolution de la tridiagonale
        solveTridiagonal(a,b,c,d);

        // Mise à jour de la grille et application des conditions aux bornes
        for(int i=1;i<M;i++) grid.set(i,j,d[i-1]);
        grid.set(0,j,lower_bc);
        grid.set(M,j,upper_bc);

        // Application de la barrière au pas de temps j
        apply_barrier_row(j);
    }

    // Retour de la valeur interpolée au strike
    int iStrike = static_cast<int>(option.getStrike()/dS);
    return grid.get(iStrike,0);
}
