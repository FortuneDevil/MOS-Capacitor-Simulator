#ifndef MYMATH_H
#define MYMATH_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "matrix.hpp"
#include "banded.hpp"
#include "block3.hpp"

// Bernoulli function B(x) = x/(e^x - 1), for any real type
template <typename T>
T bernoulli(T x){
    if (std::abs(x) < T(1e-10)) return 1 - x / 2;
    return x / std::expm1(x);
}

// Derivative of the Bernoulli function, B'(x) = B(x) * (1 - B(-x)) / x
template <typename T>
T dBernoulli(T x){
    if (std::abs(x) < T(1e-3)) return T(-0.5) + x / 6 - x * x * x / 180;
    return bernoulli(x) * (1 - bernoulli(-x)) / x;
}

double B(double x);
double dB(double x);

// Solves the tridiagonal system of a steady-state conservation law on nodes 0..N-1:
//   F(k) - F(k-1) - sink[k] * y[k] = -source[k],  with flux F(k) = up[k] * y[k+1] - down[k] * y[k]
// between nodes k and k+1, no flux on the left of node 0, and y[N-1] = last_value.
// All coefficients must be >= 0. The elimination never subtracts the (large, nearly equal) fluxes,
// so the result stays accurate even when sink and source are tiny compared to up and down.
std::vector<double> solveFluxSystem(const std::vector<double>& up, const std::vector<double>& down,
                                    const std::vector<double>& sink, const std::vector<double>& source, double last_value);

#endif
