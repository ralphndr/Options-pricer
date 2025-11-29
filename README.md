# Projet C++ : Pricing d'options financières

Ce projet implémente un système de pricing d'options financières en C++, basé sur la résolution de la PDE de Black-Scholes.  
Il couvre plusieurs types d'options et différentes méthodes numériques :

- **Options européennes** : schémas explicite, implicite et Crank-Nicolson.
- **Options américaines** : gestion de l'exercice anticipé via un solver implicite.
- Extension possible vers options exotiques et multi-actifs.

## Structure du projet

- Les fichiers `.hpp` définissent les classes et interfaces.
- Les fichiers `.cpp` implémentent la logique des classes et des solveurs.
- `main.cpp` teste l'ensemble des solveurs pour des options européennes et américaines.

## Compilation et exécution

Depuis le terminal, place-toi à la racine du dossier `Code` et utilise :

```bash
g++ src/*.cpp -Iinclude -std=c++17 -o OptionPricer.exe -mconsole
./OptionPricer.exe  
