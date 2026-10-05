# 📈 Options pricing project

This repository contains an object-oriented derivatives pricer, started as part of a programming elective at **ENSAE** and extended since. The project implements deterministic, stochastic and closed-form methods to value financial instruments, from European vanillas to exotic options, including implied-volatility skews and strike ladders under a jump-diffusion model.

---

## Project Overview

In this project, we tackle the valuation of financial derivatives within a Black-Scholes-Merton framework. The objective is to determine the fair price of an option given market parameters such as volatility, interest rates, and time to maturity.

**Key features:**
- **Dual valuation engine:** Prices options using either PDE (Partial Differential Equation) solvers or Monte Carlo simulations.
- **American option support:** Handles early exercise features via a free-boundary algorithm.
- **Exotic derivatives:** Includes pricing for Barrier, Digital, Asian, Lookback, and Chooser options.
- **Jump-diffusion smile:** Closed-form Merton pricing and a Black-Scholes implied-volatility solver to build the 1y 80–120 volatility skew, cross-checked against Monte Carlo at every strike.
- **Volatility ladder:** Marks American and exotic options at each strike's implied vol from the Merton smile, comparing flat-vol and smile-consistent prices.
- **Visualization:** Generates 3D price surfaces using Gnuplot integration.
- **High performance:** Optimized C++ architecture utilizing polymorphism and smart pointers for speed and modularity.

---

## ⚙️ How It Works

### Pricing Methods
1.  **PDE solvers:** Resolves the Black-Scholes-Merton equation on a discretized spatial-temporal grid.
    - Explicit scheme: Fast but requires respecting the CFL (stability) condition.
    - Implicit scheme: Unconditionally stable using a tridiagonal system solver.
    - Crank-Nicolson: Recommended for second-order accuracy in both time and space.
2.  **Monte Carlo simulation:** Estimates prices by averaging discounted payoffs over thousands of simulated paths.
    - PCG32 generator: Fast and reliable 32-bit random number generation.
    - Variance reduction: Implements antithetic variables to improve estimation reliability.
    - Merton Jump model: Captures market shocks by adding Poisson-distributed jumps to the price paths.
3.  **Closed-form and implied volatility:**
    - Merton (1976) price as a Poisson-weighted sum of Black-Scholes prices.
    - Black-Scholes implied volatility by bisection, inverted on out-of-the-money options.
    - Continuous-monitoring barrier formula as a reference for the PDE and Monte Carlo barrier prices.

### Market Hypotheses
The pricer assumes a frictionless market: perfect information, no arbitrage, fixed risk-free rate ($r$), and no dividends. Volatility is constant in the Black-Scholes engines; the smile comes from the Merton jump-diffusion model.

---

## 📁 Repository Structure

The project follows a modular structure to separate declarations from implementations:

```text
/Code
├── include/                # Headers (.hpp)
│   ├── PDESolver.hpp       # Abstract base class for PDE solvers
│   ├── MonteCarlo.hpp      # Stochastic simulation engine
│   ├── AmericanOption.hpp  # Specific logic for American early exercise
│   ├── Grid.hpp            # Space-time mesh management
│   └── MertonAnalytic.hpp  # Merton closed form, Black-Scholes implied vol
│
├── src/                    # Implementations (.cpp)
│   ├── main.cpp            # Entry point and coordination
│   ├── ExplicitSolver.cpp  # PDE implementation files
│   ├── MonteCarlo.cpp      # Simulation implementation files
│   └── MertonAnalytic.cpp  # Closed-form Merton and implied-vol solver
│
├── tools/                  # Analysis executables
│   ├── skew_main.cpp       # Merton implied-vol skew across strikes
│   └── ladder_main.cpp     # 1y 80-120 vol ladder across American and exotic options
│
├── output_grid.csv         # Exported price data for visualization
├── LICENSE                 # MIT License
└── README.md               # Project overview
```

## 🚀 Getting Started

It requires a minimal environment to run, with the exception of the graphical component which requires third-party software installation.

### Prerequisites
- **C++17 Compiler**: `clang++` (macOS/Linux) or `g++` version 7+.
- **Operating System**: macOS, Linux, or Windows.
- **Libraries**: Standard C++ library only.
- **Gnuplot**: Version 5.0+ is required for 3D visualization.

### Gnuplot Installation
- **macOS** (via Homebrew): `brew install gnuplot`.
- **Linux** (via APT): `sudo apt update && sudo apt install gnuplot`.
- **Windows**: Download directly from [gnuplot.info](http://www.gnuplot.info/download.html).

### Compilation
Can be compiled using a single command executed from the "Code" directory:
```bash
clang++ -std=c++17 -O2 -I./include src/*.cpp -o option_viz
```
The skew and ladder analyses are separate executables:
```bash
clang++ -std=c++17 -O2 -I./include tools/skew_main.cpp src/MertonAnalytic.cpp src/MonteCarlo.cpp -o skew
./skew --lambda 0.5 --muJ -0.1 --sigJ 0.3

clang++ -std=c++17 -O2 -I./include tools/ladder_main.cpp src/MertonAnalytic.cpp src/MonteCarlo.cpp src/CrankNicolsonSolver.cpp src/AmericanSolver.cpp src/AmericanOption.cpp src/VanillaOption.cpp src/PDESolver.cpp src/Grid.cpp -o ladder
./ladder --lambda 0.5 --muJ -0.1 --sigJ 0.3
```
## 📊 Results and Validation

The **Cpp_project_Ralph_Nader_and_Benjamin_Benisti.pdf** report details the PDE and Monte Carlo methods, with numerical tables and 3D surfaces.

### Key findings:
- **PDE solver performance**: The Crank-Nicolson scheme was found to be the most efficient, offering second-order precision while maintaining unconditional stability.
- **American vs. European**: Simulations confirmed that American puts maintain a higher value than their European counterparts due to the early exercise premium, which was calculated to be approximately 9.3% in our default scenario (6.193 vs 5.665 at S=K=100, T=1, r=5%, σ=20%).
- **Monte Carlo convergence**: While highly flexible for path-dependent options (Asian, Lookback), the Monte Carlo method was roughly 8 to 10 times slower than PDE solvers for a similar level of error.
- **Merton Jump model**: Integrating jump-diffusion processes significantly impacted the valuation (approx. 25% difference), effectively capturing market shocks that the standard Black-Scholes model misses.
- **Volatility skew**: With λ = 0.5, μ_J = −10%, σ_J = 30%, jumps steepen the 1y 80–120 skew to **2.75 vols** (ATM vol 27.7%); Monte Carlo matches the closed form within one standard error at every strike.
- **Volatility ladder**: The American put early-exercise premium rises from 3.8% (K=80) to 9.0% (K=120). A digital marked at its strike vol misses the skew contribution (0.488 flat vs 0.515 smile-consistent at the money). For an up-and-out call (B=130), closed form, PDE and Monte Carlo agree, while the jump model values it about 50% higher, since the strike vol overstates the upside diffusion. The floating-strike lookback has no strike, so it is a single mark.

---

## 👥 Authors

This project was developed by: **Ralph NADER** and **Benjamin BENISTI**.
