#include "carrier.hpp"

double carrierSteady(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p, int time, bool solve_n, bool solve_p){
    int N = V.size() - start_SC;    // semiconductor nodes, from the interface to the bulk contact
    // Coefficients of the flux between nodes i and i + 1: up * (value at i + 1) - down * (value at i)
    std::vector<double> up_n(N, 0.0), down_n(N, 0.0), up_p(N, 0.0), down_p(N, 0.0);
    // SRH recombination R = (n p - ni^2) / srh, written as sink * (own carrier) - source
    std::vector<double> sink_n(N, 0.0), sink_p(N, 0.0), source(N, 0.0);

    for (int i = 0; i < N - 1; i++){
        int gi = i + start_SC;  // global index
        double delta_x_i = x[gi + 1] - x[gi];
        double dphi = (V(gi + 1, time) - V(gi, time)) / Vt;
        double nfactor = Vt * mu_n / delta_x_i;
        double pfactor = Vt * mu_p / delta_x_i;

        up_n[i] = nfactor * B(dphi);
        down_n[i] = nfactor * B(-dphi);
        up_p[i] = pfactor * B(-dphi);
        down_p[i] = pfactor * B(dphi);

        double h = controlVolume(x, gi);
        double srh = tau_p * (n(gi, time) + ni) + tau_n * (p(gi, time) + ni);
        sink_n[i] = h * p(gi, time) / srh;
        sink_p[i] = h * n(gi, time) / srh;
        source[i] = h * ni * ni / srh;
    }

    // Solve for n, p. No flux at the interface (no current into the oxide);
    // the bulk contact keeps its equilibrium concentrations.
    // A carrier that is not solved keeps its present profile.
    std::vector<double> new_n, new_p;
    if (solve_n) new_n = solveFluxSystem(up_n, down_n, sink_n, source, n(V.size() - 1, time));
    if (solve_p) new_p = solveFluxSystem(up_p, down_p, sink_p, source, p(V.size() - 1, time));

    double max_change = 0;
    for (int i = 0; i < N - 1; i++){
        int gi = i + start_SC;  // global index
        if (solve_n){
            max_change = std::max(max_change, fabs(new_n[i] / n(gi, time) - 1));
            n(gi, time) = new_n[i];
        }
        if (solve_p){
            max_change = std::max(max_change, fabs(new_p[i] / p(gi, time) - 1));
            p(gi, time) = new_p[i];
        }
    }
    return max_change;
}
