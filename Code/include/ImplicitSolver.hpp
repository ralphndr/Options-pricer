#ifndef IMPLICIT_SOLVER_HPP
#define IMPLICIT_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

// Solveur implicite pour la résolution numérique de l'équation de Black–Scholes.
// Utilise un schéma en différences finies implicite pour garantir la stabilité
// numérique, notamment pour des pas de temps plus grands.
class ImplicitSolver : public PDESolver {
public:
    // Constructeur du solveur implicite.
    // Initialise l'option, les données de marché et les paramètres
    // de discrétisation en temps et en espace.
    ImplicitSolver(const Option& opt, const MarketData& market,
                   int M, int N, double Smax);

    // Calcule le prix de l'option par résolution implicite
    // de l'équation de Black–Scholes.
    double price() override;

private:
    // Construit les coefficients de la matrice tridiagonale
    // pour le schéma implicite à un pas de temps donné.
    void setupTridiagonal(std::vector<double>& a,
                          std::vector<double>& b,
                          std::vector<double>& c,
                          double dt);

    // Résout un système linéaire tridiagonal à l'aide de l'algorithme de Thomas.
    void solveTridiagonal(const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& c,
                          std::vector<double>& d);
};

#endif
