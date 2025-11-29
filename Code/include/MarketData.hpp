#ifndef MARKET_DATA_HPP
#define MARKET_DATA_HPP

struct MarketData {
    double r;     // Taux sans risque
    double sigma; // Volatilité
    double q;     // Dividend yield (optionnel)

    MarketData(double r, double sigma, double q = 0.0)
        : r(r), sigma(sigma), q(q) {}
};

#endif
