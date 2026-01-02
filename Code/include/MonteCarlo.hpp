#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include "MarketData.hpp"
#include "VanillaOption.hpp"
#include <string>

// Paramètres de simulation pour le Monte Carlo
struct MCParams {
    int paths = 50000;           // Nombre de trajectoires simulées
    int steps = 252;             // Nombre de pas de temps par trajectoire
    unsigned int seed = 0;       // Graine pour le générateur (0 = random_device)
    bool antithetic = true;      // Utilisation de variables antithétiques pour réduire la variance

    // Paramètres pour le modèle Merton jump-diffusion
    // Si jumpIntensity <= 0, le modèle de sauts est désactivé
    double jumpIntensity = 0.0;  // Lambda : fréquence des sauts
    double jumpMean = 0.0;       // Mu_J : moyenne logarithmique des sauts
    double jumpVol = 0.0;        // Sigma_J : volatilité des sauts

    // Paramètres pour option lookback
    std::string lookbackType = ""; // "min" ou "max", vide si non utilisé

    // Paramètres pour option chooser
    double chooseTime = -1.0;     // Fraction de T pour l'option chooser (-1 = désactivé)
};

// Résultat d'une simulation Monte Carlo
struct MCResult {
    double price;   // Prix moyen simulé
    double stderr;  // Erreur standard de l'estimation
};

// Fonction de valorisation Monte Carlo pour différentes classes d'options.
// optionClass : "Vanilla", "Digital", "Barrier", "Lookback", etc.
// optType     : Call ou Put
// spot        : prix actuel du sous-jacent
// strike      : prix d'exercice
// maturity    : maturité de l'option
// md          : données de marché
// params      : paramètres de simulation Monte Carlo
// barrierType : type de barrière pour les options barrières (facultatif)
// barrierLevel: niveau de la barrière (facultatif)
MCResult mc_price(const std::string& optionClass,
                  OptionType optType,
                  double spot,
                  double strike,
                  double maturity,
                  const MarketData& md,
                  const MCParams& params,
                  const std::string& barrierType = "",
                  double barrierLevel = -1.0);

#endif
