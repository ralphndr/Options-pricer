#include "../include/Grid.hpp"

// Constructeur de la grille pour les schémas PDE.
// Initialise les paramètres de discrétisation en espace et en temps.
// M : nombre de pas en espace (prix sous-jacent)
// N : nombre de pas en temps
// Smax : prix maximum du sous-jacent
// T : maturité de l'option
Grid::Grid(int M, int N, double Smax, double T)
    : M(M), N(N) {

    // Pas en espace et en temps
    dS = Smax / M;
    dt = T / N;

    // Initialisation de la grille avec des zéros
    // V[i][j] correspond à la valeur de l'option pour S=i*dS et t=j*dt
    V.resize(M + 1, std::vector<double>(N + 1, 0.0));
}

// Accesseur à la valeur de la grille à un point (i,j)
double Grid::get(int i, int j) const {
    return V[i][j];
}

// Modification de la valeur de la grille à un point (i,j)
void Grid::set(int i, int j, double value) {
    V[i][j] = value;
}

// Accesseurs aux tailles de la grille et aux pas
int Grid::getSizeS() const { return M; }
int Grid::getSizeT() const { return N; }
double Grid::getdS() const { return dS; }
double Grid::getdt() const { return dt; }
