#ifndef VANILLA_OPTION_HPP
#define VANILLA_OPTION_HPP

#include <algorithm>

// Type d'option : Call ou Put
enum class OptionType { Call, Put };

// Classe de base pour toutes les options.
// Définit le strike et la maturité et déclare la fonction de payoff virtuelle.
class Option {
public:
    // Constructeur de l'option
    // Initialise le prix d'exercice et la maturité
    Option(double strike, double maturity);

    // Destructeur virtuel par défaut
    virtual ~Option() = default;

    // Accesseur au prix d'exercice
    double getStrike() const;

    // Accesseur à la maturité
    double getMaturity() const;

    // Fonction de payoff purement virtuelle.
    // Chaque type d'option définit sa propre implémentation.
    virtual double payoff(double S) const = 0;

protected:
    double strike;    // Prix d'exercice
    double maturity;  // Maturité de l'option
};

// Option vanille européenne (Call ou Put)
class EuropeanOption : public Option {
public:
    // Constructeur de l'option européenne
    // Initialise le strike, la maturité et le type (Call ou Put)
    EuropeanOption(double K, double T, OptionType type);

    // Fonction de payoff pour l'option européenne
    // Call : max(S - K, 0)
    // Put  : max(K - S, 0)
    double payoff(double S) const override;

    // Indique si l'option est un Call
    bool isCall() const { return type == OptionType::Call; }

private:
    OptionType type; // Type de l'option (Call ou Put)
};

#endif
