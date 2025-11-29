#include "../include/Option.hpp"

Option::Option(double K, double T)
    : strike(K), maturity(T) {}

double Option::getStrike() const { return strike; }
double Option::getMaturity() const { return maturity; }
