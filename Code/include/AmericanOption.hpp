#ifndef AMERICAN_OPTION_HPP
#define AMERICAN_OPTION_HPP

#include "Option.hpp"

enum class AmericanOptionType { Call, Put }; // Type d’option pour AmericanOption

class AmericanOption : public Option {
public:
    AmericanOption(double strike, double maturity, AmericanOptionType type);

    double payoff(double S) const override;

private:
    AmericanOptionType optionType; // stocke le type de l'option
};

#endif
