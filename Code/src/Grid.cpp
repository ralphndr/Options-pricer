#include "../include/Grid.hpp"

Grid::Grid(int M, int N, double Smax, double T)
    : M(M), N(N), Smax(Smax), T(T) {

    dS = Smax / M;
    dt = T / N;

    V.resize(M + 1, std::vector<double>(N + 1, 0.0));
}

double Grid::get(int i, int j) const {
    return V[i][j];
}

void Grid::set(int i, int j, double value) {
    V[i][j] = value;
}

int Grid::getSizeS() const { return M; }
int Grid::getSizeT() const { return N; }
double Grid::getdS() const { return dS; }
double Grid::getdt() const { return dt; }
