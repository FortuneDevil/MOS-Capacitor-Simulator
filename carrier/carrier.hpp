#ifndef CARRIER_H
#define CARRIER_H
#include "../main/main.hpp"

// Steady-state carrier continuity (Scharfetter-Gummel + SRH recombination) for column `time` of V, n, p.
// Each returns the largest relative change made to n or p.
double carrierSteady(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int time);
double carrierEq(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p);
double carrierDC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p);
double carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond);

// Small-signal AC analysis around the DC solution (column 1 of V, n, p) at angular frequency omega.
// Solves the linearised Poisson and continuity equations together for the phasors V_ac, n_ac, p_ac
// produced by a gate signal of 1 V amplitude. Returns the complex capacitance per unit area.
std::complex<double> smallSignalAC(const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, double omega,
                                   std::vector<std::complex<double>>& V_ac, std::vector<std::complex<double>>& n_ac,
                                   std::vector<std::complex<double>>& p_ac);

#endif
