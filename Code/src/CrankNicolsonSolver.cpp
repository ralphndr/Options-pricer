#include "../include/CrankNicolsonSolver.hpp"
#include "../include/VanillaOption.hpp"
#include "../include/ExoticOption.hpp"
#include <vector>
#include <iostream>

// Constructeur du solveur Crank–Nicolson.
// Initialise l'option, les données de marché et la grille via le constructeur PDESolver.
CrankNicolsonSolver::CrankNicolsonSolver(const Option& opt,
                                         const MarketData& market,
                                         int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

// Calcule le prix de l'option via le schéma de Crank–Nicolson
// pour différents types d'options (vanilla, digitale, barrière).
double CrankNicolsonSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dt = grid.getdt();
    double dS = grid.getdS();
    double r = data.r;

    // Initialisation de la grille à maturité avec le payoff de l'option
    for(int i=0;i<=M;++i)
        grid.set(i,N,option.payoff(i*dS));

    std::vector<double> a(M-1), b(M-1), c(M-1), d(M-1);
    double T = N * dt;  // Temps total jusqu'à l'échéance

    // Identification du type d'option pour les conditions aux bords
    bool is_call = true;
    if (const EuropeanOption* euro_opt = dynamic_cast<const EuropeanOption*>(&option)) {
        is_call = euro_opt->isCall();
    } else if (const BarrierOption* bo = dynamic_cast<const BarrierOption*>(&option)) {
        is_call = (bo->getOptionType() == OptionType::Call);
    }

    // Boucle arrière en temps
    for(int j=N-1;j>=0;j--){
        double tau = j * dt;  
        double time_to_expiry = T - tau;
        
        // Détermination des conditions aux bords pour ce pas de temps
        double lower_bc, upper_bc;
        double S_max = M * dS;

        if (const DigitalOption* dig = dynamic_cast<const DigitalOption*>(&option)) {
            // Options digitales : payoff fixe avec actualisation
            bool call = (dig->getOptionType() == OptionType::Call);
            if (call) {
                lower_bc = 0.0;
                upper_bc = dig->getPayout() * std::exp(-r * time_to_expiry);
            } else {
                lower_bc = dig->getPayout() * std::exp(-r * time_to_expiry);
                upper_bc = 0.0;
            }
        } else {
            // Options vanilla et autres : conditions asymptotiques
            if(is_call) {
                lower_bc = 0.0;
                upper_bc = S_max - option.getStrike() * std::exp(-r * time_to_expiry);
            } else {
                lower_bc = option.getStrike() * std::exp(-r * time_to_expiry);
                upper_bc = 0.0;
            }
        }
        
        // Construction des coefficients CN pour les points intérieurs
        for(int i=1;i<M;i++){
            double S = i*dS;
            double alpha = 0.25*dt*(data.sigma*data.sigma*S*S/dS/dS - data.r*S/dS);
            double beta  = -0.5*dt*(data.sigma*data.sigma*S*S/dS/dS + data.r);
            double gamma = 0.25*dt*(data.sigma*data.sigma*S*S/dS/dS + data.r*S/dS);

            a[i-1] = -alpha;
            b[i-1] = 1 - beta;
            c[i-1] = -gamma;

            // Côté droit (explicit) pour Crank–Nicolson
            d[i-1] = alpha*grid.get(i-1,j+1) + (1+beta)*grid.get(i,j+1) + gamma*grid.get(i+1,j+1);
            
            // Ajustement des RHS pour les conditions aux bords
            if(i==1) d[i-1] -= a[i-1] * lower_bc;
            if(i==M-1) d[i-1] -= c[i-1] * upper_bc;
        }

        // Résolution tridiagonale par l'algorithme de Thomas
        std::vector<double> c_star(M-1), d_star(M-1);
        c_star[0] = c[0]/b[0];
        d_star[0] = d[0]/b[0];
        for(int i=1;i<M-1;i++){
            double m = b[i] - a[i]*c_star[i-1];
            c_star[i] = c[i]/m;
            d_star[i] = (d[i] - a[i]*d_star[i-1])/m;
        }

        d[M-2] = d_star[M-2];
        for(int i=M-3;i>=0;i--){
            d[i] = d_star[i] - c_star[i]*d[i+1];
        }

        // Mise à jour de la grille et application des conditions aux bords
        for(int i=1;i<M;i++) grid.set(i,j,d[i-1]);
        grid.set(0,j,lower_bc);
        grid.set(M,j,upper_bc);
    }

    // Post-traitement : application des conditions de barrière après résolution PDE
    if (const BarrierOption* bo = dynamic_cast<const BarrierOption*>(&option)) {
        bool knockIn = bo->isKnockIn();
        bool call = (bo->getOptionType() == OptionType::Call);
        double B = bo->getBarrierLevel();
        
        // Boucle sur tous les points de la grille pour appliquer la barrière
        for (int j = N; j >= 0; --j) {
            for (int i = 0; i <= M; ++i) {
                double S = i * dS;
                bool hitBarrier;
                
                if (call) {
                    // Call : barrière en haut
                    hitBarrier = (S >= B);
                } else {
                    // Put : barrière en bas
                    hitBarrier = (S <= B);
                }
                
                if (knockIn) {
                    // Knock-in : valeur nulle si barrière NON atteinte
                    if (!hitBarrier) grid.set(i, j, 0.0);
                } else {
                    // Knock-out : valeur nulle si barrière atteinte
                    if (hitBarrier) grid.set(i, j, 0.0);
                }
            }
        }
    }

    // Retour de la valeur interpolée au strike
    int iStrike = static_cast<int>(option.getStrike()/dS);
    return grid.get(iStrike,0);
}
