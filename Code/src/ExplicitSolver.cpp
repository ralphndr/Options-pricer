#include "../include/ExplicitSolver.hpp"
#include "../include/VanillaOption.hpp"
#include "../include/ExoticOption.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>

// Constructeur du solveur explicite.
// Initialise l'option, les données de marché et la grille via le constructeur PDESolver.
ExplicitSolver::ExplicitSolver(const Option& opt,
                               const MarketData& market,
                               int M, int N, double Smax)
    : PDESolver(opt, market, M, N, Smax) {}

// Initialisation de la grille à maturité et application des conditions de barrière.
// Cette fonction prépare le payoff et applique les conditions aux bords.
void ExplicitSolver::initializeBoundary() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();

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

    // Initialisation à maturité (payoff)
    for(int i = 0; i <= M; ++i) {
        double S = i * dS;
        grid.set(i, N, option.payoff(S));
    }

    // Application de la barrière à maturité
    apply_barrier_row(N);

    // Conditions aux bornes à maturité
    grid.set(0, N, option.payoff(0.0));
    grid.set(M, N, option.payoff(M * dS));
}

// Calcule le prix de l'option via le schéma explicite.
double ExplicitSolver::price() {
    int M = grid.getSizeS();
    int N = grid.getSizeT();
    double dS = grid.getdS();
    double dt = grid.getdt();
    double r = data.r;
    double sigma = data.sigma;
    double T = N * dt;  // Temps total jusqu'à l'échéance
    double Smax = M * dS;

    // Vérification de la condition CFL
    // lambda(S) = 0.5 * sigma^2 * S^2 * dt / (dS^2) <= 0.5
    // Le point critique est S = Smax
    double lambda_max = 0.5 * sigma * sigma * Smax * Smax * dt / (dS * dS);
    if (lambda_max > 0.5) {
        std::cerr << "Erreur CFL: lambda_max = " << lambda_max << " > 0.5\n";
        std::cerr << "Besoin N >= " << static_cast<int>(std::ceil(T * sigma * sigma * M * M / 0.5)) << "\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Préparation de la grille et application des conditions aux bords
    initializeBoundary();

    // Détermination du type d'option pour les conditions aux bords
    bool is_call = true;
    if (const EuropeanOption* euro_opt = dynamic_cast<const EuropeanOption*>(&option)) {
        is_call = euro_opt->isCall();
    } else if (const BarrierOption* bo = dynamic_cast<const BarrierOption*>(&option)) {
        is_call = (bo->getOptionType() == OptionType::Call);
    }

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

    // Boucle arrière en temps pour la résolution explicite
    for(int j = N - 1; j >= 0; --j) {
        double tau = j * dt;  // Temps écoulé
        double time_to_expiry = T - tau;

        // Calcul des points intérieurs de la grille
        for(int i = 1; i < M; ++i) {
            double S = i * dS;
            double V_right = grid.get(i+1, j+1);
            double V_center = grid.get(i, j+1);
            double V_left = grid.get(i-1, j+1);
            
            // Dérivées par différences centrées
            double dV_dS = (V_right - V_left) / (2.0 * dS);
            double d2V_dS2 = (V_right - 2.0*V_center + V_left) / (dS * dS);
            
            // PDE explicite: V^n_i = V^{n+1}_i + dt * [0.5*σ²*S²*∂²V/∂S² + r*S*∂V/∂S - r*V]
            double V = V_center + dt * (0.5 * sigma*sigma * S*S * d2V_dS2 
                                       + r * S * dV_dS 
                                       - r * V_center);
            grid.set(i, j, V);
        }

        // Conditions aux bornes en fonction du type d'option
        double S_max = M * dS;
        if (const DigitalOption* dig = dynamic_cast<const DigitalOption*>(&option)) {
            bool call = (dig->getOptionType() == OptionType::Call);
            if (call) {
                grid.set(0, j, 0.0);
                grid.set(M, j, dig->getPayout() * std::exp(-r * time_to_expiry));
            } else {
                grid.set(0, j, dig->getPayout() * std::exp(-r * time_to_expiry));
                grid.set(M, j, 0.0);
            }
        } else {
            if(is_call) {
                grid.set(0,j, 0.0);
                grid.set(M,j, S_max - option.getStrike() * std::exp(-r * time_to_expiry));
            } else {
                grid.set(0,j, option.getStrike() * std::exp(-r * time_to_expiry));
                grid.set(M,j, 0.0);
            }
        }

        // Application de la barrière à ce pas de temps
        apply_barrier_row(j);
    }

    // Retour de la valeur interpolée au strike
    int iStrike = static_cast<int>(option.getStrike() / dS);
    return grid.get(iStrike, 0);
}
