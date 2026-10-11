#include "carrier.hpp"

// Equilibrium (VG = 0) is the steady-state problem solved on column 0
double carrierEq(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p){
    return carrierSteady(x, V, n, p, 0);
}
