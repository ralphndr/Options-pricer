#ifndef EXOTIC_OPTION_HPP
#define EXOTIC_OPTION_HPP

#include "VanillaOption.hpp"
#include <string>

// Classe représentant une option barrière (knock-in ou knock-out).
// Le type de barrière (up/down, in/out) est spécifié sous forme de chaîne,
// et le payoff est évalué comme celui d’une option européenne conditionnellement
// à l’activation ou non de la barrière, gérée au niveau du solveur.
class BarrierOption : public Option {
public:
    // Constructeur de l’option barrière.
    // Initialise le strike, la maturité, le type de l’option,
    // le type de barrière et son niveau.
    BarrierOption(double K, double T, OptionType t, const std::string& bType, double bLevel)
        : Option(K, T), type(t), barrierType(bType), barrierLevel(bLevel) {}

    // Fonction de payoff à maturité.
    // Le payoff correspond à celui d’une option européenne vanille,
    // sous réserve que la condition de barrière soit satisfaite.
    double payoff(double S) const override {
        if (type == OptionType::Call) return std::max(0.0, S - strike);
        return std::max(0.0, strike - S);
    }

    // Accesseur au type de l’option (Call ou Put)
    OptionType getOptionType() const { return type; }

    // Accesseur au type de barrière
    std::string getBarrierType() const { return barrierType; }

    // Accesseur au niveau de la barrière
    double getBarrierLevel() const { return barrierLevel; }

    // Indique si l’option est de type knock-in
    bool isKnockIn() const {
        return barrierType == "downin" || barrierType == "upin" || barrierType == "knockin";
    }

    // Indique si l’option est de type knock-out
    bool isKnockOut() const {
        return barrierType == "downout" || barrierType == "upout" || barrierType == "knockout";
    }

private:
    // Type de l’option (Call ou Put)
    OptionType type;

    // Type de barrière (up/down, in/out)
    std::string barrierType;

    // Niveau de la barrière
    double barrierLevel;
};

// Classe représentant une option digitale (cash-or-nothing).
// Le payoff est discontinu et correspond à un montant fixe
// si l’option est dans la monnaie à maturité.
class DigitalOption : public Option {
public:
    // Constructeur de l’option digitale.
    // Initialise le strike, la maturité, le type de l’option
    // et le montant du paiement fixe.
    DigitalOption(double K, double T, OptionType t, double cash)
        : Option(K, T), type(t), payout(cash) {}

    // Fonction de payoff à maturité.
    // Le paiement est égal au montant fixe si la condition
    // d’exercice est satisfaite, nul sinon.
    double payoff(double S) const override {
        if (type == OptionType::Call) return (S > strike) ? payout : 0.0;
        return (S < strike) ? payout : 0.0;
    }

    // Accesseur au type de l’option
    OptionType getOptionType() const { return type; }

    // Accesseur au montant du paiement fixe
    double getPayout() const { return payout; }

private:
    // Type de l’option (Call ou Put)
    OptionType type;

    // Montant du paiement fixe
    double payout;
};

#endif
