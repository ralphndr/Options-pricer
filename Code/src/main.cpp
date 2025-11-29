#include <iostream>
#include "MarketData.hpp"
#include "EuropeanOption.hpp"
#include "ExplicitSolver.hpp"
#include "ImplicitSolver.hpp"
#include "CrankNicolsonSolver.hpp"
#include "AmericanOption.hpp"
#include "AmericanSolver.hpp"

int main() {
    MarketData market(0.05, 0.2);
    EuropeanOption call(100.0, 1.0, OptionType::Call);

    ExplicitSolver explicitSolver(call, market, 200, 5000, 150.0);
    std::cout << "Prix call (explicite) : " << explicitSolver.price() << std::endl;

    ImplicitSolver implicitSolver(call, market, 200, 5000, 150.0);
    std::cout << "Prix call (implicite) : " << implicitSolver.price() << std::endl;

    CrankNicolsonSolver cnSolver(call, market, 200, 5000, 150.0);
    std::cout << "Prix call (Crank-Nicolson) : " << cnSolver.price() << std::endl;

    // Test du solver américain
    AmericanOption amCall(100.0, 1.0, AmericanOptionType::Call);
    AmericanSolver amSolver(amCall, market, 200, 5000, 150.0);
    std::cout << "Prix call américain : " << amSolver.price() << std::endl;

    return 0;
}




