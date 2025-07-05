#ifndef CARRIER_H
#define CARRIER_H
#include "../main/main.hpp"

void carrierEq(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int startSC);
void carrierDC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int startSC);
void carrierAC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int startSC);
void carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond, int startSC);

#endif