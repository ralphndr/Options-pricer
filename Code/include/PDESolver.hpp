#ifndef PDE_SOLVER_HPP
#define PDE_SOLVER_HPP

#include "VanillaOption.hpp"
#include "MarketData.hpp"
#include "Grid.hpp"

// Classe de base pour les solveurs PDE (différences finies) 
// pour l'équation de Black–Scholes.
// Fournit l'accès à l'option, aux données de marché et à la grille de calcul.
class PDESolver {
public:
    // Constructeur du solveur PDE.
    // Initialise l'option, les données de marché et la grille
    // en fonction du nombre de points en espace et en temps et de Smax.
    PDESolver(const Option& opt, const MarketData& market,
              int M, int N, double Smax);

    // Destructeur virtuel par défaut
    virtual ~PDESolver() = default;

    // Méthode virtuelle pure pour calculer le prix de l'option
    virtual double price() = 0;

    // Accès en lecture à la grille calculée
    const Grid& getGrid() const;

    // Accès en écriture à la grille calculée
    Grid& getGrid();

protected:
    const Option& option;   // Option à valoriser
    const MarketData& data; // Données de marché
    Grid grid;              // Grille pour les différences finies
};

#endif
