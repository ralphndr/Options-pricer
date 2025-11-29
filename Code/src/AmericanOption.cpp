#include "../include/AmericanOption.hpp"
#include <algorithm>

AmericanOption::AmericanOption(double strike, double maturity, AmericanOptionType type)
    : Option(strike, maturity), optionType(type) {} // constructeur Option valide

double AmericanOption::payoff(double S) const {
    if(optionType == AmericanOptionType::Call)
        return std::max(S - strike, 0.0);
    else
        return std::max(strike - S, 0.0);
}
