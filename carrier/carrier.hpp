#ifndef CARRIER_H
#define CARRIER_H
#include "../main/main.hpp"

void carrierEq(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p);
void carrierDC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p);
void carrierAC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p);
void carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond);

#endif