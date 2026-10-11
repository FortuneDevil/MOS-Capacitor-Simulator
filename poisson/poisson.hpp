#ifndef POISSON_H
#define POISSON_H

#include "../main/main.hpp"

// One Newton step of Poisson's equation. n and p follow V exponentially (fixed quasi-Fermi levels).
// Returns the largest change made to V.
double poissonSolver(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond);

#endif
