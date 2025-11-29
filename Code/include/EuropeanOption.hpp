#ifndef EUROPEAN_OPTION_HPP
#define EUROPEAN_OPTION_HPP

#include "Option.hpp"
#include <algorithm>

enum class OptionType { Call, Put };

class EuropeanOption : public Option {
public:
    EuropeanOption(double K, double T, OptionType type);

    double payoff(double S) const override;

private:
    OptionType type;
};

#endif
