#ifndef CARRIER_H
#define CARRIER_H
#include "../main/main.hpp"

// Steady-state carrier continuity (Scharfetter-Gummel + SRH recombination) for column `time` of V, n, p.
// A carrier that is not solved keeps its present profile. Returns the largest relative change made to n or p.
double carrierSteady(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int time,
                     bool solve_n = true, bool solve_p = true);

// Carrier solve for a steady-state condition (EQUILIBRIUM, DC or DEEP_DEPLETION)
double carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond);

#endif
