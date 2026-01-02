#ifndef MARKET_DATA_HPP
#define MARKET_DATA_HPP

// Structure contenant les données de marché nécessaires à la valorisation
// d'une option via les équations de Black–Scholes.
struct MarketData {
    double r;     // Taux d'intérêt sans risque
    double sigma; // Volatilité du sous-jacent
    double q;     // Dividend yield (facultatif, par défaut 0)

    // Constructeur de la structure MarketData.
    // Initialise les paramètres de marché pour le calcul du prix de l'option.
    MarketData(double r, double sigma, double q = 0.0)
        : r(r), sigma(sigma), q(q) {}
};

#endif
