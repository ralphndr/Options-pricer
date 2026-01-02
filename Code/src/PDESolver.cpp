#include "../include/PDESolver.hpp"

// PDESolver : classe de base pour résoudre les PDE d'options
// Variables :
//   option : référence constante vers l'option à valoriser
//   data   : référence constante vers les données de marché (taux, volatilité)
//   grid   : grille S x t pour stocker les valeurs de l'option
PDESolver::PDESolver(const Option& opt,
                     const MarketData& market,
                     int M, int N, double Smax)
    : option(opt), data(market), grid(M, N, Smax, opt.getMaturity()) {}

// Accesseur constant vers la grille (lecture seule)
const Grid& PDESolver::getGrid() const {
    return grid;
}

// Accesseur non-constant vers la grille (lecture/écriture)
Grid& PDESolver::getGrid() {
    return grid;
}
