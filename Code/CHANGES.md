# Changes after the course submission

## Fixes
1. **Crank-Nicolson barrier options** (`src/CrankNicolsonSolver.cpp`): the knock-out
   condition is now imposed at every time step instead of once after the full
   backward solve. Before: CN up-and-out call (K=100, B=130) = 10.29, identical to the
   vanilla. After: 3.35 (implicit 3.36; continuous-monitoring closed form 3.33;
   MC 3.41 with 2000 steps, 3.56 with 252 steps due to discrete monitoring).
2. **American solver boundary terms** (`src/AmericanSolver.cpp`): the boundary values
   at the new time level are now added to the right-hand side of the implicit system.
   Before: American call 10.2056 < European call 10.2917 (impossible). After: equal.
3. **Monte Carlo standard error with antithetic variates** (`src/MonteCarlo.cpp`):
   computed from pair averages, since the two paths of a pair are not independent.
   ATM call, 100k paths: SE 0.0468 without antithetics, 0.0328 with them.
4. Added missing `#include <cmath>` so the code also compiles with g++ on Linux.

## New: Merton implied-vol skew (`include/MertonAnalytic.hpp`, `src/MertonAnalytic.cpp`, `tools/skew_main.cpp`)
Closed-form Merton price, Black-Scholes implied vol by bisection, 80-120 strike ladder
at T=1, each strike cross-checked against the Monte Carlo engine with jumps.
With the report's jump parameters (lambda=0.5, muJ=-0.10, sigJ=0.30): ATM vol 27.7%,
80-120 skew 2.75 vols, MC within 1 standard error of the closed form at every strike.

    clang++ -std=c++17 -O2 -I./include tools/skew_main.cpp src/MertonAnalytic.cpp src/MonteCarlo.cpp -o skew
    ./skew --lambda 0.5 --muJ -0.1 --sigJ 0.3

## Note on evaluation point
PDE solvers return the grid value at index int(K/dS). Choose Smax so K is a grid node
(e.g. Smax=200 with even M); otherwise the price is read slightly below S=K.

## New: 1y 80-120 vol ladder across American and exotic options (`tools/ladder_main.cpp`)
For each strike, read the implied vol off the Merton smile and mark every product
struck there at that vol: European and American put (PDE, early-exercise premium),
digital call (flat vol vs smile-consistent -dC/dK), up-and-out call B=130 (closed form,
PDE and MC as independent checks, plus the Merton MC price), Asian call (MC at the
strike vol vs Merton). The floating-strike lookback has no strike, so it is one mark.

    clang++ -std=c++17 -O2 -I./include tools/ladder_main.cpp src/MertonAnalytic.cpp src/MonteCarlo.cpp src/CrankNicolsonSolver.cpp src/AmericanSolver.cpp src/AmericanOption.cpp src/VanillaOption.cpp src/PDESolver.cpp src/Grid.cpp -o ladder
    ./ladder --lambda 0.5 --muJ -0.1 --sigJ 0.3
