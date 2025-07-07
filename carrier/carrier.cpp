#include "carrier.hpp"

void carrier(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, Condition cond){
    switch (cond) {
        case EQUILIBRIUM:
            carrierEq(x, V, n, p);
            break;
        case DC:
            carrierDC(x, V, n, p);
            break;
        case AC:
            carrierAC(x, V, n, p);
            break;
        default:
            break;
    }
}