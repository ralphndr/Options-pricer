#ifndef AMERICAN_OPTION_HPP
#define AMERICAN_OPTION_HPP

#include "VanillaOption.hpp"

// Type d'option américaine : Call ou Put
enum class AmericanOptionType { Call, Put };

// Classe représentant une option américaine vanille.
// Elle définit la structure du produit et le calcul du payoff.
// La gestion de l'exercice anticipé est déléguée aux solveurs numériques.
class AmericanOption : public Option {
public:
    // Constructeur de la classe AmericanOption.
    // Initialise le prix d'exercice, la maturité et le type de l'option.
    AmericanOption(double strike, double maturity, AmericanOptionType type);

    // Fonction de payoff de l'option.
    // Calcule la valeur intrinsèque de l'option en fonction du prix du sous-jacent.
    double payoff(double S) const override;

    // Accesseur au type de l'option (Call ou Put)
    AmericanOptionType getType() const { return optionType; }

private:
    // Type de l'option américaine
    AmericanOptionType optionType;
};

#endif
