#ifndef AC_H
#define AC_H
#include "../main/main.hpp"

// Small-signal AC analysis around the DC solution in column `time` of V, n, p, at angular frequency omega.
// Solves the linearised Poisson and continuity equations together for the phasors V_ac, n_ac, p_ac
// produced by a gate signal of 1 V amplitude. Returns the complex capacitance per unit area C_complex:
// the admittance is j omega C_complex, so C = Re(C_complex) and G/omega = -Im(C_complex).
std::complex<double> smallSignalAC(const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, int time, double omega,
                                   std::vector<std::complex<double>>& V_ac, std::vector<std::complex<double>>& n_ac,
                                   std::vector<std::complex<double>>& p_ac);

#endif
