#include "../include/MertonAnalytic.hpp"
#include <cmath>
#include <algorithm>

static double ncdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double bs_price_eu(double S, double K, double T, double r, double q,
                   double sigma, OptionType type) {
    double dfR = std::exp(-r * T), dfQ = std::exp(-q * T);
    if (T <= 0.0 || sigma <= 0.0) {               // degenerate: intrinsic on forward
        double fwdIntr = S * dfQ - K * dfR;
        return (type == OptionType::Call) ? std::max(0.0, fwdIntr) : std::max(0.0, -fwdIntr);
    }
    double sd = sigma * std::sqrt(T);
    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / sd;
    double d2 = d1 - sd;
    if (type == OptionType::Call) return S * dfQ * ncdf(d1) - K * dfR * ncdf(d2);
    return K * dfR * ncdf(-d2) - S * dfQ * ncdf(-d1);
}

double merton_price(double S, double K, double T, double r, double q,
                    double sigma, double lambda, double muJ, double sigJ,
                    OptionType type) {
    // kappa = E[J - 1]: expected relative jump size
    double kappa   = std::exp(muJ + 0.5 * sigJ * sigJ) - 1.0;
    double lamP    = lambda * (1.0 + kappa);          // lambda' (jump-adjusted intensity)
    double logW    = -lamP * T;                       // log Poisson weight for n = 0
    double price   = 0.0;

    for (int n = 0; n < 200; ++n) {
        if (n > 0) logW += std::log(lamP * T) - std::log(static_cast<double>(n));
        double w = std::exp(logW);
        // Conditional on n jumps: lognormal with extra variance n*sigJ^2,
        // and drift shifted so the discounted price stays a martingale.
        double sigN = std::sqrt(sigma * sigma + n * sigJ * sigJ / T);
        double rN   = r - lambda * kappa + n * std::log(1.0 + kappa) / T;
        // Hull's form: with lambda' weights, plain BS at rate rN is exact
        // (equivalent to lambda weights times exp((rN - r) T) * BS(rN)).
        price += w * bs_price_eu(S, K, T, rN, q, sigN, type);
        if (n > lamP * T && w < 1e-14) break;         // past the Poisson mode and negligible
    }
    return price;
}

double bs_implied_vol(double price, double S, double K, double T,
                      double r, double q, OptionType type) {
    // No-arbitrage bounds: price must lie between intrinsic value and the upper bound
    double lo = bs_price_eu(S, K, T, r, q, 1e-8, type);
    double hi = (type == OptionType::Call) ? S * std::exp(-q * T) : K * std::exp(-r * T);
    if (price < lo - 1e-12 || price > hi) return -1.0;

    double a = 1e-6, b = 5.0;                         // vol bracket: 0.0001% to 500%
    for (int it = 0; it < 200; ++it) {                // BS price is increasing in vol (vega > 0)
        double m = 0.5 * (a + b);
        if (bs_price_eu(S, K, T, r, q, m, type) < price) a = m; else b = m;
        if (b - a < 1e-10) break;
    }
    return 0.5 * (a + b);
}
