#include "../include/MonteCarlo.hpp"
#include <random>
#include <cmath>
#include <algorithm>

// CDF de la loi normale standard
static double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// Prix Black-Scholes pour call/put européen (utilisé au choix du chooser)
static double bs_price(double S, double K, double tau, double r, double q, double sigma, OptionType type) {
    if (tau <= 0.0) {
        return (type == OptionType::Call) ? std::max(0.0, S - K) : std::max(0.0, K - S);
    }
    double volSqrt = sigma * std::sqrt(tau);
    if (volSqrt <= 0.0) {
        return (type == OptionType::Call) ? std::max(0.0, S * std::exp(-q * tau) - K * std::exp(-r * tau))
                                          : std::max(0.0, K * std::exp(-r * tau) - S * std::exp(-q * tau));
    }
    double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * tau) / volSqrt;
    double d2 = d1 - volSqrt;
    double df_r = std::exp(-r * tau);
    double df_q = std::exp(-q * tau);
    double call = S * df_q * norm_cdf(d1) - K * df_r * norm_cdf(d2);
    if (type == OptionType::Call) return call;
    return call - S * df_q + K * df_r; // put-call parity
}

// Implémentation minimale du générateur PCG32 pour Monte Carlo
// Inspiré par O'Neill, générateur pseudo-aléatoire 32 bits de haute qualité
class PCG32 {
public:
    PCG32(uint64_t graine, uint64_t seq = 1) { reseed(graine, seq); }

    void reseed(uint64_t graine, uint64_t seq) {
        state = 0U;
        inc = (seq << 1u) | 1u;
        next();
        state += graine;
        next();
    }

    uint32_t next() {
        uint64_t etatAncien = state;
        state = etatAncien * 6364136223846793005ULL + inc;
        uint32_t xorshifted = static_cast<uint32_t>(((etatAncien >> 18u) ^ etatAncien) >> 27u);
        uint32_t rot = static_cast<uint32_t>(etatAncien >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    // Génère un nombre uniforme dans (0,1)
    double uniforme() {
        uint32_t u = next();
        const double denom = 1.0 / 4294967296.0; // 2^32
        double x = (static_cast<double>(u) + 0.5) * denom;
        if (x <= 0.0) x = 1e-12;
        if (x >= 1.0) x = 1.0 - 1e-12;
        return x;
    }

    // Génère un nombre suivant la loi normale standard
    double normale() {
        if (hasSpare) {
            hasSpare = false;
            return spare;
        }
        double u1 = uniforme();
        double u2 = uniforme();
        double r = std::sqrt(-2.0 * std::log(u1));
        double theta = 2.0 * M_PI * u2;
        spare = r * std::sin(theta);
        hasSpare = true;
        return r * std::cos(theta);
    }

private:
    uint64_t state = 0;
    uint64_t inc = 0;
    double spare = 0.0;
    bool hasSpare = false;
};

// Vérifie si la barrière est franchie selon le type de barrière
static bool barriere_franchie(const std::string& type, double S, double B) {
    if (type == "downout" || type == "downin") return S <= B;
    if (type == "upout" || type == "upin") return S >= B;
    return false;
}

// Calcul du prix d'une option par simulation Monte Carlo
// Gère : Européenne, Asiatique, Lookback, Chooser, Barrière
MCResult mc_price(const std::string& classeOption,
                  OptionType typeOpt,
                  double spot,
                  double strike,
                  double maturite,
                  const MarketData& md,
                  const MCParams& params,
                  const std::string& typeBarriere,
                  double niveauBarriere) {

    const int nbChemins = std::max(1, params.paths); // nombre de chemins simulés
    const int nbEtapes = std::max(1, params.steps);  // nombre de pas de temps

    // Initialisation du générateur aléatoire
    uint64_t graine = params.seed;
    if (graine == 0) {
        std::random_device rd;
        graine = (static_cast<uint64_t>(rd()) << 32) ^ rd();
    }
    PCG32 rng(graine, 1u);

    const double dt = maturite / static_cast<double>(nbEtapes);  // pas de temps
    const double vol = md.sigma * std::sqrt(dt);                 // volatilité ajustée

    // Paramètres pour Jump-Diffusion (Merton)
    const bool utiliserSauts = (params.jumpIntensity > 0.0 && params.jumpVol >= 0.0);
    const double lambda = utiliserSauts ? params.jumpIntensity : 0.0;
    const double kappa = utiliserSauts ? (std::exp(params.jumpMean + 0.5 * params.jumpVol * params.jumpVol) - 1.0) : 0.0;
    const double drift = (md.r - md.q - 0.5 * md.sigma * md.sigma - lambda * kappa) * dt;
    const double facteurActualisation = std::exp(-md.r * maturite);

    double somme = 0.0;       // somme des payoffs actualisés
    double sommeCarres = 0.0; // somme des carrés pour erreur type

    // Fonction payoff pour option européenne
    auto payoff_europeen = [&](double S) {
        return (typeOpt == OptionType::Call) ? std::max(0.0, S - strike)
                                             : std::max(0.0, strike - S);
    };

    // Générateur poisson pour sauts
    std::mt19937 rngSaut(static_cast<uint32_t>(graine ^ 0x9e3779b97f4a7c15ULL));
    std::poisson_distribution<int> pois(utiliserSauts ? lambda * dt : 0.0);

    // simuler un chemin avec tableau de nombres aléatoires générés
    auto simulate_path = [&](const std::vector<std::vector<double>>& zs, 
                              const std::vector<std::vector<int>>& nbSautsVec,
                              const std::vector<std::vector<double>>& zSautsVec,
                              int idxChemin, bool useAnti) -> double {
        double S_courant = spot;
        double sommeS = 0.0;
        double S_min = spot;
        double S_max = spot;
        bool barriereFranchieFlag = false;
        OptionType typeChoisi = typeOpt;
        int etapeChoisie = -1;
        if (params.chooseTime > 0.0)
            etapeChoisie = static_cast<int>(params.chooseTime * nbEtapes);

        for (int etape = 0; etape < nbEtapes; ++etape) {
            double z_base = zs[etape][idxChemin];
            double z_utilise = useAnti ? -z_base : z_base;

            double chocSaut = 0.0;
            if (utiliserSauts) {
                int nbSauts = nbSautsVec[etape][idxChemin];
                if (nbSauts > 0) {
                    double z_saut_base = zSautsVec[etape][idxChemin];
                    double z_saut = useAnti ? -z_saut_base : z_saut_base;
                    chocSaut = nbSauts * params.jumpMean + std::sqrt(static_cast<double>(nbSauts)) * params.jumpVol * z_saut;
                }
            }

            S_courant *= std::exp(drift + vol * z_utilise + chocSaut);

            if (!typeBarriere.empty() && barriere_franchie(typeBarriere, S_courant, niveauBarriere))
                barriereFranchieFlag = true;

            if (classeOption == "asian") sommeS += S_courant;

            if (!params.lookbackType.empty()) {
                S_min = std::min(S_min, S_courant);
                S_max = std::max(S_max, S_courant);
            }

            if (etape == etapeChoisie && params.chooseTime > 0.0) {
                double tauRestant = maturite - (etape + 1) * dt;
                double valeurCall = bs_price(S_courant, strike, std::max(tauRestant, 0.0), md.r, md.q, md.sigma, OptionType::Call);
                double valeurPut  = bs_price(S_courant, strike, std::max(tauRestant, 0.0), md.r, md.q, md.sigma, OptionType::Put);
                typeChoisi = (valeurCall >= valeurPut) ? OptionType::Call : OptionType::Put;
            }
        }

        double payoffFinal = 0.0;
        if (!params.lookbackType.empty()) {
            if (params.lookbackType == "min") payoffFinal = std::max(0.0, S_courant - S_min);
            else payoffFinal = std::max(0.0, S_max - S_courant);
        } else if (params.chooseTime > 0.0) {
            payoffFinal = (typeChoisi == OptionType::Call) ? std::max(0.0, S_courant - strike)
                                                            : std::max(0.0, strike - S_courant);
        } else if (classeOption == "asian") {
            double moyenne = sommeS / static_cast<double>(nbEtapes);
            payoffFinal = (typeOpt == OptionType::Call) ? std::max(0.0, moyenne - strike)
                                                        : std::max(0.0, strike - moyenne);
        } else {
            payoffFinal = payoff_europeen(S_courant);
        }

        if (!typeBarriere.empty()) {
            bool knockIn = (typeBarriere.find("in") != std::string::npos);
            if (knockIn && !barriereFranchieFlag) payoffFinal = 0.0;
            if (!knockIn && barriereFranchieFlag) payoffFinal = 0.0;
        }

        return facteurActualisation * payoffFinal;
    };

    // Pré-génération des nombres aléatoires pour toutes les étapes
    // Pour antithetic: générer nbChemins/2 nombres de base, réutiliser 2 fois (z et -z)
    // Pour normal: générer nbChemins nombres
    int nbZsToGenerate = params.antithetic ? (nbChemins + 1) / 2 : nbChemins;
    std::vector<std::vector<double>> zs(nbEtapes, std::vector<double>(nbZsToGenerate));
    std::vector<std::vector<int>> nbSautsVec(nbEtapes, std::vector<int>(nbZsToGenerate, 0));
    std::vector<std::vector<double>> zSautsVec(nbEtapes, std::vector<double>(nbZsToGenerate, 0.0));

    for (int etape = 0; etape < nbEtapes; ++etape) {
        for (int i = 0; i < nbZsToGenerate; ++i) {
            zs[etape][i] = rng.normale();
            if (utiliserSauts) {
                nbSautsVec[etape][i] = pois(rngSaut);
                if (nbSautsVec[etape][i] > 0) {
                    zSautsVec[etape][i] = rng.normale();
                }
            }
        }
    }

    // Simulation des chemins
    if (params.antithetic) {
        // Antithétique: chaque paire réutilise les mêmes nombres z et -z
        int nbPaires = (nbChemins + 1) / 2;
        for (int p = 0; p < nbPaires; ++p) {
            // Chemin 1: utilise z directement
            double payoff1 = simulate_path(zs, nbSautsVec, zSautsVec, p, false);
            somme += payoff1;
            sommeCarres += payoff1 * payoff1;

            // Chemin 2: utilise -z (réutilise le même "p" pour les nombres, mais avec signe opposé)
            if (2 * p + 1 < nbChemins) {
                double payoff2 = simulate_path(zs, nbSautsVec, zSautsVec, p, true);
                somme += payoff2;
                sommeCarres += payoff2 * payoff2;
            }
        }
    } else {
        // Sans antithétique: chemins indépendants, chacun utilise un ensemble de nombres unique
        for (int i = 0; i < nbChemins; ++i) {
            double payoffActualise = simulate_path(zs, nbSautsVec, zSautsVec, i, false);
            somme += payoffActualise;
            sommeCarres += payoffActualise * payoffActualise;
        }
    }

    // Moyenne et erreur type
    double moyenne = somme / static_cast<double>(nbChemins);
    double variance = (sommeCarres / nbChemins) - moyenne * moyenne;
    if (variance < 0.0) variance = 0.0;
    double erreurType = std::sqrt(variance / nbChemins);

    return { moyenne, erreurType };
}
