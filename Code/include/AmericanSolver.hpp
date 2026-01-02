#ifndef AMERICANSOLVER_HPP
#define AMERICANSOLVER_HPP

#include "AmericanOption.hpp"
#include "PDESolver.hpp"

// Solveur numérique pour la valorisation d'une option américaine.
// La spécificité américaine (exercice anticipé) est prise en compte
// via une condition de type obstacle lors de la résolution de l'équation
// de Black–Scholes.
class AmericanSolver : public PDESolver {
public:
    // Constructeur du solveur pour option américaine.
    // Initialise l'option, les données de marché et les paramètres
    // de discrétisation en temps et en espace.
    AmericanSolver(const AmericanOption& opt, const MarketData& market,
                   int M, int N, double Smax);

    // Calcule le prix de l'option américaine par résolution numérique
    // de l'équation de Black–Scholes avec condition d'exercice anticipé.
    double price() override;
};

#endif
