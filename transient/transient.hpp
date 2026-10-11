#ifndef TRANSIENT_H
#define TRANSIENT_H

#include "../main/main.hpp"

// Time-dependent simulation. The state lives in the TRANSIENT column of V, n, p, together with the
// quasi-Fermi potentials Phi (column 0: electrons, column 1: holes), defined by
//   n = ni exp((V - Phi_n) / Vt),   p = ni exp((Phi_p - V) / Vt)

// Sets Phi from the V, n, p that are in the TRANSIENT column (the starting state)
void transientInit(const Matrix& V, const Matrix& n, const Matrix& p, Matrix& Phi);

// Advances the state by one time step dt (backward Euler) with the gate node held at V_gate.
// Poisson and both continuity equations are solved together by Newton's method.
// Returns the number of Newton iterations, or -1 if it did not converge.
int transientStep(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Matrix& Phi, double V_gate, double dt);

#endif
