#include "carrier.hpp"

double carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond){
    switch (cond) {
        case EQUILIBRIUM:
            return carrierEq(x, V, n, p);
        case DC:
            return carrierDC(x, V, n, p);
        default:
            return 0;
    }
}
