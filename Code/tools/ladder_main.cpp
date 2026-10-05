// 1y 80-120 vol ladder across American and exotic options.
//
// Idea ("marking" by strike): the Merton model gives a smile. For each strike K
// we read the implied vol sigma_K off that smile, then mark every product struck
// at K with sigma_K using our Black-Scholes engines (PDE and MC). Where possible
// we compare against a smile-consistent or model price to show what a flat
// per-strike vol misses.
//
// Build (from Code/):
//   clang++ -std=c++17 -O2 -I./include tools/ladder_main.cpp src/MertonAnalytic.cpp
//     src/MonteCarlo.cpp src/CrankNicolsonSolver.cpp src/AmericanSolver.cpp
//     src/AmericanOption.cpp src/VanillaOption.cpp src/PDESolver.cpp src/Grid.cpp -o ladder
//   (one line)
// Run:
//   ./ladder --lambda 0.5 --muJ -0.1 --sigJ 0.3

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "MertonAnalytic.hpp"
#include "MonteCarlo.hpp"
#include "MarketData.hpp"
#include "VanillaOption.hpp"
#include "ExoticOption.hpp"
#include "AmericanOption.hpp"
#include "AmericanSolver.hpp"
#include "CrankNicolsonSolver.hpp"

// Closed-form continuously monitored up-and-out call, valid for B > K
// (Reiner-Rubinstein / Haug: A - B + C - D). Used as the reference for the PDE and MC.
static double uo_call_exact(double S, double K, double B, double T, double r, double sig) {
    auto N = [](double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); };
    double sq = sig * std::sqrt(T), mu = (r - 0.5 * sig * sig) / (sig * sig), d = std::exp(-r * T);
    double x  = std::log(S / K) / sq + (1 + mu) * sq,  x1 = std::log(S / B) / sq + (1 + mu) * sq;
    double y  = std::log(B * B / (S * K)) / sq + (1 + mu) * sq, y1 = std::log(B / S) / sq + (1 + mu) * sq;
    double A  = S * N(x)  - K * d * N(x - sq);
    double Bt = S * N(x1) - K * d * N(x1 - sq);
    double C  = S * std::pow(B / S, 2 * (mu + 1)) * N(-y)  - K * d * std::pow(B / S, 2 * mu) * N(-y + sq);
    double D  = S * std::pow(B / S, 2 * (mu + 1)) * N(-y1) - K * d * std::pow(B / S, 2 * mu) * N(-y1 + sq);
    return A - Bt + C - D;
}

// PDE solvers fill the whole grid; read the value at the spot (not at the strike).
// Smax and M are chosen so the spot sits exactly on a grid node.
static double value_at_spot(PDESolver& solver, double spot) {
    solver.price();
    const Grid& g = solver.getGrid();
    int i = static_cast<int>(std::lround(spot / g.getdS()));
    return g.get(i, 0);
}

int main(int argc, char** argv) {
    double S = 100.0, T = 1.0, r = 0.05, q = 0.0, sigma = 0.20;
    double lambda = 0.5, muJ = -0.10, sigJ = 0.30;    // report's jump parameters
    double B = 130.0;                                  // up-and-out barrier
    int paths = 100000, steps = 252;
    unsigned seed = 42;
    for (int i = 1; i + 1 < argc; i += 2) {
        if      (!std::strcmp(argv[i], "--lambda")) lambda = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--muJ"))    muJ    = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--sigJ"))   sigJ   = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--sigma"))  sigma  = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--barrier"))B      = std::atof(argv[i + 1]);
        else if (!std::strcmp(argv[i], "--paths"))  paths  = std::atoi(argv[i + 1]);
    }
    const int M = 800, N = 2000;
    const double Smax = 400.0;                         // dS = 0.5 -> spot 100 is node 200
    const double dt = T / steps;
    const double bgk = 0.5826;                         // Broadie-Glasserman-Kou constant

    std::printf("1y vol ladder | S=%.0f r=%.2f | Merton: sigma=%.2f lambda=%.2f muJ=%.2f sigJ=%.2f | UO barrier=%.0f\n\n",
                S, r, sigma, lambda, muJ, sigJ, B);
    std::printf("%5s %7s | %9s %9s %6s | %9s %9s | %9s %9s %9s %9s | %9s %9s\n",
                "K", "vol(%)", "Euro put", "Amer put", "EEP%", "Dig flat", "Dig smile",
                "UO exact", "UO PDE", "UO MC-BS", "UO MC-Mrt", "Asn MC-BS", "Asn Mrt");

    for (double K = 80.0; K <= 120.0 + 1e-9; K += 5.0) {
        // 1) Strike vol from the Merton smile (OTM option for the inversion)
        OptionType otm = (K < S) ? OptionType::Put : OptionType::Call;
        double volK = bs_implied_vol(merton_price(S, K, T, r, q, sigma, lambda, muJ, sigJ, otm),
                                     S, K, T, r, q, otm);
        MarketData mdK(r, volK, q);

        // 2) European vs American put marked at vol_K (PDE)
        EuropeanOption euroPut(K, T, OptionType::Put);
        CrankNicolsonSolver cnEuro(euroPut, mdK, M, N, Smax);
        double euro = value_at_spot(cnEuro, S);
        AmericanOption amPut(K, T, AmericanOptionType::Put);
        AmericanSolver amSolver(amPut, mdK, M, N, Smax);
        double amer = value_at_spot(amSolver, S);

        // 3) Digital call: flat vol_K (PDE) vs smile-consistent price.
        //    A digital is minus the strike-derivative of the call price, so the
        //    smile-consistent value is -dC/dK of the Merton call (central difference).
        DigitalOption dig(K, T, OptionType::Call, 1.0);
        CrankNicolsonSolver cnDig(dig, mdK, M, N, Smax);
        double digFlat = value_at_spot(cnDig, S);
        const double h = 0.01;
        double digSmile = -(merton_price(S, K + h, T, r, q, sigma, lambda, muJ, sigJ, OptionType::Call)
                          - merton_price(S, K - h, T, r, q, sigma, lambda, muJ, sigJ, OptionType::Call)) / (2 * h);

        // 4) Up-and-out call: PDE at vol_K vs MC at vol_K (independent check) vs MC under Merton.
        //    MC monitors the barrier daily; shifting it down by exp(-0.5826 sigma sqrt(dt))
        //    makes the discrete MC comparable to the (near-)continuous PDE barrier.
        BarrierOption uo(K, T, OptionType::Call, "upout", B);
        CrankNicolsonSolver cnUO(uo, mdK, M, N, Smax);
        double uoPDE = value_at_spot(cnUO, S);

        MCParams mpBS; mpBS.paths = paths; mpBS.steps = steps; mpBS.seed = seed;
        double uoMcBS = mc_price("european", OptionType::Call, S, K, T, mdK, mpBS, "upout",
                                 B * std::exp(-bgk * volK * std::sqrt(dt))).price;

        MarketData mdDiff(r, sigma, q);               // Merton: diffusion vol + jumps
        MCParams mpJ = mpBS; mpJ.jumpIntensity = lambda; mpJ.jumpMean = muJ; mpJ.jumpVol = sigJ;
        double uoMcMrt = mc_price("european", OptionType::Call, S, K, T, mdDiff, mpJ, "upout",
                                  B * std::exp(-bgk * sigma * std::sqrt(dt))).price;

        // 5) Asian (arithmetic, daily fixings) call: MC at vol_K vs MC under Merton
        double asnBS  = mc_price("asian", OptionType::Call, S, K, T, mdK,    mpBS).price;
        double asnMrt = mc_price("asian", OptionType::Call, S, K, T, mdDiff, mpJ).price;

        std::printf("%5.0f %7.2f | %9.4f %9.4f %6.2f | %9.4f %9.4f | %9.4f %9.4f %9.4f %9.4f | %9.4f %9.4f\n",
                    K, 100 * volK, euro, amer, 100 * (amer / euro - 1), digFlat, digSmile,
                    uo_call_exact(S, K, B, T, r, volK), uoPDE, uoMcBS, uoMcMrt, asnBS, asnMrt);
    }

    // Floating-strike lookback has no strike, so it is a single mark, not a ladder.
    MCParams mpL; mpL.paths = paths; mpL.steps = steps; mpL.seed = seed; mpL.lookbackType = "min";
    double atmVol = bs_implied_vol(merton_price(S, S, T, r, q, sigma, lambda, muJ, sigJ, OptionType::Call),
                                   S, S, T, r, q, OptionType::Call);
    double lbBS = mc_price("european", OptionType::Call, S, S, T, MarketData(r, atmVol, q), mpL).price;
    MCParams mpLJ = mpL; mpLJ.jumpIntensity = lambda; mpLJ.jumpMean = muJ; mpLJ.jumpVol = sigJ;
    double lbMrt = mc_price("european", OptionType::Call, S, S, T, MarketData(r, sigma, q), mpLJ).price;
    std::printf("\nFloating-strike lookback call (no strike): MC at ATM vol %.2f%% = %.4f | MC under Merton = %.4f\n",
                100 * atmVol, lbBS, lbMrt);
    std::printf("\nColumns: EEP = early-exercise premium of the American put over the European, both at vol_K.\n");
    std::printf("Dig flat = digital marked at the strike's vol; Dig smile = -dC/dK under Merton (includes skew).\n");
    std::printf("UO exact = continuous-monitoring closed form at vol_K; UO PDE and UO MC-BS are independent checks\n");
    std::printf("(PDE monitors at its %d time steps, so it sits slightly above; daily MC with the BGK shift sits close).\n", N);
    std::printf("UO MC-Mrt = same option under the jump model: lower diffusion vol, mostly downward jumps.\n");
    return 0;
}
