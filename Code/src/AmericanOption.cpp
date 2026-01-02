#include "../include/AmericanOption.hpp"
#include <algorithm> // pour std::max

// Constructeur de la classe AmericanOption.
// Initialise les paramètres communs via le constructeur de la classe de base Option
// (prix d'exercice et maturité), puis définit le type de l'option (Call ou Put).
AmericanOption::AmericanOption(double strike, double maturity, AmericanOptionType type)
    : Option(strike, maturity), optionType(type) {}

// Fonction de payoff pour une option américaine.
// Le payoff est identique à celui d'une option européenne correspondante :
// - Call : max(S - K, 0)
// - Put  : max(K - S, 0)
// La spécificité américaine (exercice anticipé) est traitée
// au niveau du solveur numérique, et non dans cette fonction locale.
double AmericanOption::payoff(double S) const {
    if(optionType == AmericanOptionType::Call)
        return std::max(S - strike, 0.0);
    else
        return std::max(strike - S, 0.0);
}
