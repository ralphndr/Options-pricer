#include "../include/EuropeanOption.hpp"

EuropeanOption::EuropeanOption(double K, double T, OptionType type)
    : Option(K, T), type(type) {}

double EuropeanOption::payoff(double S) const {
    if (type == OptionType::Call)
        return std::max(S - strike, 0.0);
    else
        return std::max(strike - S, 0.0);
}
