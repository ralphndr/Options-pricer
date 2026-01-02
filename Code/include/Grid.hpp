#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

// Classe représentant une grille discrète pour les différences finies.
// Permet de stocker et manipuler les valeurs de l'option en fonction
// du prix du sous-jacent et du temps.
class Grid {
public:
    // Constructeur de la grille.
    // Initialise une grille de taille M x N pour l'espace et le temps.
    // M : nombre de points en espace (prix du sous-jacent)
    // N : nombre de points en temps
    // Smax : valeur maximale du sous-jacent
    // T : maturité de l'option
    Grid(int M, int N, double Smax, double T);

    // Retourne la valeur stockée à la position (i,j) sur la grille
    double get(int i, int j) const;

    // Modifie la valeur stockée à la position (i,j) sur la grille
    void set(int i, int j, double value);

    // Retourne le nombre de points en espace
    int getSizeS() const;

    // Retourne le nombre de points en temps
    int getSizeT() const;

    // Retourne le pas en espace (prix du sous-jacent)
    double getdS() const;

    // Retourne le pas en temps
    double getdt() const;

private:
    int M, N;                            // Dimensions de la grille (espace et temps)
    double dS, dt;                       // Pas en espace et en temps
    std::vector<std::vector<double>> V;  // Matrice des valeurs de l'option
};

#endif
