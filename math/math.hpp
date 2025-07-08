#ifndef MYMATH_H
#define MYMATH_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "matrix.hpp"

double B(double x);

Matrix solver(Matrix &A, Matrix &b);

double max_relative_error(const Matrix& A, const Matrix& B, int col);
#endif