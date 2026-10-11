#include "poisson.hpp"
#include <cstdlib>

double poissonSolver(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond){
    int N = V.size(); //  V is a column matrix with V.rows = N
    int time = (cond == DC) ? 1 : 0;
    BandMatrix<double> J(N, 1, 1);
    std::vector<double> F(N, 0.0);

    // Gate (V = VG, already stored in V) and bulk (V = 0) are fixed
    J(0, 0) = 1.0;
    J(N - 1, N - 1) = 1.0;

    for (int i = 1; i < N - 1; i++){
        double delta_x_minus = x[i] - x[i - 1];
        double delta_x_i = x[i + 1] - x[i];
        double e_minus = (i <= start_SC) ? e_ox : e_si; // permittivity between nodes i-1 and i
        double e_plus = (i < start_SC) ? e_ox : e_si;   // permittivity between nodes i and i+1
        double h = controlVolume(x, i);                 // zero in the oxide, so no charge there

        J(i, i - 1) = e_minus / delta_x_minus;
        J(i, i + 1) = e_plus / delta_x_i;
        J(i, i) = -(e_minus / delta_x_minus + e_plus / delta_x_i) - q * h * (n(i, time) + p(i, time)) / Vt;

        F[i] = -(e_minus * (V(i - 1, time) - V(i, time)) / delta_x_minus
               + e_plus * (V(i + 1, time) - V(i, time)) / delta_x_i
               + q * h * (ND - NA + p(i, time) - n(i, time)));
    }

    std::vector<double> delta_v = solveBanded(J, F);

    double max_change = 0;
    for (int i = 1; i < N - 1; i++){
        double dv = delta_v[i];
        if (std::isnan(dv)){
            std::cerr << "delta_v[" << i << "] is NaN\n";
            exit(1);
        }
        // Damp large steps so that the exponentials in n, p don't overshoot
        if (fabs(dv) > Vt) dv = (dv > 0 ? 1 : -1) * Vt * log(1 + fabs(dv) / Vt);

        V(i, time) += dv;
        n(i, time) *= exp(dv / Vt);
        p(i, time) *= exp(-dv / Vt);
        max_change = std::max(max_change, fabs(dv));
    }
    return max_change;
}
