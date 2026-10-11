#include "carrier.hpp"

double carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond){
    bool p_type = NA > ND;
    switch (cond) {
        case EQUILIBRIUM:
        case DC:
            return carrierSteady(x, V, n, p, column(cond));
        case DEEP_DEPLETION:
            // Only the majority carrier reaches steady state, the minority carrier stays frozen
            return carrierSteady(x, V, n, p, column(cond), !p_type, p_type);
        default:
            return 0;
    }
}
