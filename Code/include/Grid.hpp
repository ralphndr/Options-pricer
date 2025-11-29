#ifndef GRID_HPP
#define GRID_HPP

#include <vector>

class Grid {
public:
    Grid(int M, int N, double Smax, double T);

    double get(int i, int j) const;
    void set(int i, int j, double value);

    int getSizeS() const;
    int getSizeT() const;
    double getdS() const;
    double getdt() const;

private:
    int M, N;
    double Smax, T;
    double dS, dt;
    std::vector<std::vector<double>> V;
};

#endif
