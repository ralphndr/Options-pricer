# 📈 Options pricing project

This repository contains an object-oriented derivatives pricer developed as part of a programming elective at **ENSAE**. The project implements both deterministic and stochastic methods to value financial instruments, ranging from European vanillas to complex exotic options.

---

## Project Overview

In this project, we tackle the valuation of financial derivatives within a Black-Scholes-Merton framework. The objective is to determine the fair price of an option given market parameters such as volatility, interest rates, and time to maturity.

**Key features:**
- **Dual valuation engine:** Prices options using either PDE (Partial Differential Equation) solvers or Monte Carlo simulations.
- **American option support:** Handles early exercise features via a free-boundary algorithm.
- **Exotic derivatives:** Includes pricing for Barrier, Digital, Asian, Lookback, and Chooser options.
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

### Market Hypotheses
The pricer assumes a frictionless market: perfect information, no arbitrage, constant volatility ($\sigma$), fixed risk-free rate ($r$), and no dividends.

---

## 📁 Repository Structure

The project follows a modular structure to separate declarations from implementations:

```text
/Code
├── include/                # Headers (.hpp)
│   ├── PDESolver.hpp       # Abstract base class for PDE solvers
│   ├── MonteCarlo.hpp      # Stochastic simulation engine
│   ├── AmericanOption.hpp  # Specific logic for American early exercise
│   └── Grid.hpp            # Space-time mesh management
│
├── src/                    # Implementations (.cpp)
│   ├── main.cpp            # Entry point and coordination
│   ├── ExplicitSolver.cpp  # PDE implementation files
│   └── MonteCarlo.cpp      # Simulation implementation files
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
## 📊 Results and Validation

The project includes a comprehensive analysis of the performance and accuracy of each pricing method. For a complete view of the numerical tables and 3D surfaces, please refer to the **Cpp_project_Ralph_Nader_and_Benjamin_Benisti.pdf** file included in this repository.

### Key findings:
- **PDE solver performance**: The Crank-Nicolson scheme was found to be the most efficient, offering second-order precision while maintaining unconditional stability.
- **American vs. European**: Simulations confirmed that American puts maintain a higher value than their European counterparts due to the early exercise premium, which was calculated to be approximately 9.4% in our default scenario.
- **Monte Carlo convergence**: While highly flexible for path-dependent options (Asian, Lookback), the Monte Carlo method was roughly 8 to 10 times slower than PDE solvers for a similar level of error.
- **Merton Jump model**: Integrating jump-diffusion processes significantly impacted the valuation (approx. 25% difference), effectively capturing market shocks that the standard Black-Scholes model misses.

---

## 👥 Authors

This project was developed by: **Ralph NADER** and **Benjamin BENISTI**.
