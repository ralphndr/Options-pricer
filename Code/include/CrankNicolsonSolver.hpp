#ifndef CRANK_NICOLSON_SOLVER_HPP
#define CRANK_NICOLSON_SOLVER_HPP

#include "PDESolver.hpp"
#include <vector>

// Solveur numérique basé sur le schéma de Crank–Nicolson
// pour la résolution de l'équation de Black–Scholes.
class CrankNicolsonSolver : public PDESolver {
public:
    // Constructeur du solveur Crank–Nicolson.
    // Initialise l'option, les données de marché et les paramètres
    // de discrétisation en temps et en espace.
    CrankNicolsonSolver(const Option& opt, const MarketData& market,
                        int M, int N, double Smax);

    // Calcule le prix de l'option par résolution numérique
    // de l'équation de Black–Scholes à l'aide du schéma de Crank–Nicolson.
    double price() override;

private:
    // Construit les coefficients de la matrice tridiagonale
    // associés au schéma de Crank–Nicolson pour un pas de temps donné.
    void setupTridiagonalCN(std::vector<double>& a,
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
