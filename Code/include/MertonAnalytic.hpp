#ifndef MERTON_ANALYTIC_HPP
#define MERTON_ANALYTIC_HPP

#include "VanillaOption.hpp"   // OptionType

// Black-Scholes price of a European call/put (continuous dividend yield q).
double bs_price_eu(double S, double K, double T, double r, double q,
                   double sigma, OptionType type);

// Merton (1976) jump-diffusion price, closed form.
// Log-jump sizes ~ N(muJ, sigJ^2), jump arrivals ~ Poisson(lambda * T).
// Price = sum over n jumps of Poisson weight * Black-Scholes price with
// jump-adjusted volatility and drift.
double merton_price(double S, double K, double T, double r, double q,
                    double sigma, double lambda, double muJ, double sigJ,
                    OptionType type);

// Black-Scholes implied volatility by bisection.
// Returns -1 if the price lies outside the no-arbitrage bounds.
double bs_implied_vol(double price, double S, double K, double T,
                      double r, double q, OptionType type);

#endif
