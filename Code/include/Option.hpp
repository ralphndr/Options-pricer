#ifndef OPTION_HPP
#define OPTION_HPP

class Option {
public:
    Option(double strike, double maturity);
    virtual ~Option() = default;

    double getStrike() const;
    double getMaturity() const;

    // Payoff virtuel → chaque type d'option définit le sien
    virtual double payoff(double S) const = 0;

protected:
    double strike;
    double maturity;
};

#endif
