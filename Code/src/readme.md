# Sources du projet

Ce dossier contient tous les fichiers sources (`.cpp`) qui implémentent la logique des classes et solveurs définis dans `include/` :

- `main.cpp` : programme principal, crée les options, instancie les solveurs et affiche les prix.  
- `MarketData.cpp` : implémentation de la classe `MarketData`.  
- `EuropeanOption.cpp` et `AmericanOption.cpp` : implémentation des méthodes de payoff.  
- `ExplicitSolver.cpp`, `ImplicitSolver.cpp`, `CrankNicolsonSolver.cpp`, `AmericanSolver.cpp` : implémentation des solveurs.  
- `Grid.cpp` : implémentation de la grille de discrétisation.

## Compilation et exécution

Depuis la racine `/Code` :

```bash
g++ src/*.cpp -Iinclude -std=c++17 -o OptionPricer.exe -mconsole
./OptionPricer.exe

