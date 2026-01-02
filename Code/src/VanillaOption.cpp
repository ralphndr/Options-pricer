#include "../include/VanillaOption.hpp"

// Variables :
//   strike   : prix d'exercice
//   maturity : maturité de l'option (en années)
Option::Option(double K, double T)
    : strike(K), maturity(T) {}

double Option::getStrike() const { return strike; }
double Option::getMaturity() const { return maturity; }

// Variables :
//   type : Call ou Put
EuropeanOption::EuropeanOption(double K, double T, OptionType type)
    : Option(K, T), type(type) {}

// Payoff européen classique
double EuropeanOption::payoff(double S) const {
    return (type == OptionType::Call) ? std::max(S - strike, 0.0)
                                     : std::max(strike - S, 0.0);
}
