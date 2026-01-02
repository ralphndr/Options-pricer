#ifndef EXPLICIT_SOLVER_HPP
#define EXPLICIT_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

// Solveur explicite pour la résolution numérique de l'équation de Black–Scholes.
// Utilise un schéma en différences finies explicite.
// Les conditions aux bords et le payoff sont pris en compte dans la discrétisation.
class ExplicitSolver : public PDESolver {
public:
    // Constructeur du solveur explicite.
    // Initialise l'option, les données de marché et les paramètres
    // de discrétisation en temps et en espace.
    ExplicitSolver(const Option& opt, const MarketData& market,
                   int M, int N, double Smax);

    // Calcule le prix de l'option par résolution explicite
    // de l'équation de Black–Scholes.
    double price() override;

private:
    // Initialise les conditions aux bords pour la grille des prix
    // du sous-jacent, nécessaires à l'évolution explicite.
    void initializeBoundary();
};

#endif
