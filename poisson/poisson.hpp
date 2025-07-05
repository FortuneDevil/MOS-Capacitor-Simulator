#ifndef POISSON_H
#define POISSON_H

#include "../main/main.hpp"

void poissonSolver(const std::vector<double>& x, Matrix& V, const Matrix& n, const Matrix& p, Condition cond);

#endif