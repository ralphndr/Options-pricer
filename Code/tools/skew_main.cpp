// Merton jump-diffusion implied-vol skew across a 1y strike ladder.
// For each strike: closed-form Merton price -> BS implied vol,
// cross-checked against the Monte Carlo engine (mc_price with jumps).
//
// Build (from Code/):
//   clang++ -std=c++17 -O2 -I./include tools/skew_main.cpp src/MertonAnalytic.cpp src/MonteCarlo.cpp -o skew
// Run:
//   ./skew                         (default parameters)
//   ./skew --lambda 1 --muJ -0.15  (override any parameter)

// System headers first: on macOS <cstdio> defines `stderr` as a macro, and
// MCResult has a member called `stderr`, so the macro must be visible
// consistently (same include order as main.cpp).
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "MertonAnalytic.hpp"
#include "MonteCarlo.hpp"
#include "MarketData.hpp"

int main(int argc, char** argv) {
    // Market and model parameters (chosen before running, not fitted)
    double S = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    double lambda = 0.5, muJ = -0.10, sigJ = 0.10;   // ~1 jump every 2 years, avg -10%
    int paths = 400000;
    unsigned seed = 42;

    for (int i = 1; i + 1 < argc; i += 2) {
        if      (!std::strcmp(argv[i], "--spot"))   S      = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--T"))      T      = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--r"))      r      = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--sigma"))  sigma  = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--lambda")) lambda = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--muJ"))    muJ    = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--sigJ"))   sigJ   = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--paths"))  paths  = std::atoi(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--seed"))   seed   = static_cast<unsigned>(std::atoi(argv[i + 1]));
    }

    MarketData md(r, sigma, q);
    MCParams mp;
    mp.paths = paths;
    mp.steps = 12;          // European payoff: log-increments are exact, step count doesn't bias
    mp.seed = seed;
    mp.antithetic = true;
    mp.jumpIntensity = lambda;
    mp.jumpMean = muJ;
    mp.jumpVol = sigJ;

    std::printf("Merton skew | S=%.0f T=%.2f r=%.3f sigma=%.3f lambda=%.2f muJ=%.3f sigJ=%.3f\n",
                S, T, r, sigma, lambda, muJ, sigJ);
    std::printf("%6s %4s %12s %10s %12s %10s %9s %10s\n",
                "K", "type", "Merton", "IV(%)", "MC", "MC se", "z-score", "MC IV(%)");

    std::vector<double> strikes;
    for (double K = 80.0; K <= 120.0 + 1e-9; K += 5.0) strikes.push_back(K);

    double iv80 = 0.0, iv100 = 0.0, iv120 = 0.0;
    for (double K : strikes) {
        // Out-of-the-money option at each strike: puts below spot, calls above.
        // Same implied vol either way (put-call parity), but OTM prices are
        // pure time value, so the vol inversion is better conditioned.
        OptionType type = (K < S) ? OptionType::Put : OptionType::Call;

        double px  = merton_price(S, K, T, r, q, sigma, lambda, muJ, sigJ, type);
        double iv  = bs_implied_vol(px, S, K, T, r, q, type);

        MCResult mc = mc_price("european", type, S, K, T, md, mp);
        double ivMc = bs_implied_vol(mc.price, S, K, T, r, q, type);
        double z    = (mc.price - px) / mc.stderr;

        std::printf("%6.0f %4s %12.5f %10.3f %12.5f %10.5f %9.2f %10.3f\n",
                    K, type == OptionType::Put ? "P" : "C",
                    px, 100 * iv, mc.price, mc.stderr, z, 100 * ivMc);

        if (std::fabs(K - 80.0)  < 1e-9) iv80  = iv;
        if (std::fabs(K - 100.0) < 1e-9) iv100 = iv;
        if (std::fabs(K - 120.0) < 1e-9) iv120 = iv;
    }

    std::printf("\nATM vol: %.3f%%  |  80-120 skew (IV80 - IV120): %.2f vols  |  diffusion-only vol: %.1f%%\n",
                100 * iv100, 100 * (iv80 - iv120), 100 * sigma);
    return 0;
}
