#include "transient.hpp"

constexpr int max_newton = 100;
constexpr double tol_newton = 1e-10;    // Largest change of V, Phi_n, Phi_p at convergence (V)

void transientInit(const Matrix& V, const Matrix& n, const Matrix& p, Matrix& Phi){
    int time = column(TRANSIENT);
    for (int i = start_SC; i < V.size(); i++){
        Phi(i, 0) = V(i, time) - Vt * log(n(i, time) / ni);
        Phi(i, 1) = V(i, time) + Vt * log(p(i, time) / ni);
    }
}

int transientStep(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Matrix& Phi, double V_gate, double dt){
    int time = column(TRANSIENT);
    int N = V.size() - start_SC;    // semiconductor nodes, from the interface to the bulk contact
    double C_ox = e_ox / (x[start_SC] - x[0]);

    // Concentrations at the previous time
    std::vector<double> n_old(N), p_old(N);
    for (int i = 0; i < N; i++){
        n_old[i] = n(i + start_SC, time);
        p_old[i] = p(i + start_SC, time);
    }

    int iter = 0;
    bool converged = false;
    std::vector<double> n_new(N), p_new(N);
    while (iter < max_newton && !converged){
        for (int i = 0; i < N; i++){
            int gi = i + start_SC;  // global index
            n_new[i] = ni * exp((V(gi, time) - Phi(gi, 0)) / Vt);
            p_new[i] = ni * exp((Phi(gi, 1) - V(gi, time)) / Vt);
        }

        // The unknowns at node i are X_i = (V, Phi_n, Phi_p). The equations of node i are
        //   F(i) - F(i-1) + S_i = 0
        // with F(i) = (displacement, -electron flux, -hole flux) between nodes i and i + 1 and S_i the
        // local terms. Newton's method solves the linearised system
        //   up dX_(i+1) - down dX_i - (up dX_i - down dX_(i-1))[element i-1] - local dX_i = -residual
        std::vector<Mat3<double>> up(N, zeroMat3<double>()), down(N, zeroMat3<double>()), local(N, zeroMat3<double>());
        std::vector<Vec3<double>> flux(N, zeroVec3<double>()), residual(N, zeroVec3<double>());

        for (int i = 0; i < N - 1; i++){
            int gi = i + start_SC;  // global index
            double delta_x_i = x[gi + 1] - x[gi];
            double dphi = (V(gi + 1, time) - V(gi, time)) / Vt;
            double B_plus = B(dphi), B_minus = B(-dphi);
            double ratio = dB(dphi) / B_plus;   // B'(dphi) / B(dphi)

            // Poisson: e_si * (V_(i+1) - V_i) / delta_x
            flux[i][0] = e_si * (V(gi + 1, time) - V(gi, time)) / delta_x_i;
            up[i][0][0] = down[i][0][0] = e_si / delta_x_i;

            // Electron flux nfactor * [B(dphi) n_(i+1) - B(-dphi) n_i], written with the difference of
            // quasi-Fermi potentials so that drift and diffusion are never subtracted from each other
            double nfactor = Vt * mu_n / delta_x_i;
            double d_n = (Phi(gi + 1, 0) - Phi(gi, 0)) / Vt;
            double F_n = (d_n > 0) ? nfactor * B_minus * n_new[i] * expm1(-d_n)
                                   : - nfactor * B_plus * n_new[i + 1] * expm1(d_n);
            flux[i][1] = - F_n;
            up[i][1][0] = - F_n * (1 + ratio) / Vt;
            up[i][1][1] = nfactor * B_plus * n_new[i + 1] / Vt;
            down[i][1][0] = - F_n * ratio / Vt;
            down[i][1][1] = nfactor * B_minus * n_new[i] / Vt;

            // Hole flux pfactor * [B(dphi) p_i - B(-dphi) p_(i+1)], written the same way
            double pfactor = Vt * mu_p / delta_x_i;
            double d_p = (Phi(gi + 1, 1) - Phi(gi, 1)) / Vt;
            double F_p = (d_p > 0) ? pfactor * B_minus * p_new[i + 1] * expm1(-d_p)
                                   : - pfactor * B_plus * p_new[i] * expm1(d_p);
            flux[i][2] = - F_p;
            up[i][2][0] = - F_p * ratio / Vt;
            up[i][2][2] = pfactor * B_minus * p_new[i + 1] / Vt;
            down[i][2][0] = - F_p * (1 + ratio) / Vt;
            down[i][2][2] = pfactor * B_plus * p_new[i] / Vt;
        }

        for (int i = 0; i < N - 1; i++){
            int gi = i + start_SC;  // global index
            double h = controlVolume(x, gi);
            double n0 = n_new[i], p0 = p_new[i];
            double dn = h * n0 / Vt, dp = h * p0 / Vt;  // change of h * n (or h * p) per volt

            // SRH recombination. n p - ni^2 = ni^2 (exp((Phi_p - Phi_n) / Vt) - 1)
            double srh = tau_p * (n0 + ni) + tau_n * (p0 + ni);
            double R = ni * ni * expm1((Phi(gi, 1) - Phi(gi, 0)) / Vt) / srh;
            double r_n = (p0 - R * tau_p) / srh;
            double r_p = (n0 - R * tau_n) / srh;

            residual[i][0] = flux[i][0] + q * h * (ND - NA + p0 - n0);
            residual[i][1] = flux[i][1] + h * (R + (n0 - n_old[i]) / dt);
            residual[i][2] = flux[i][2] - h * (R + (p0 - p_old[i]) / dt);
            if (i > 0){
                for (int k = 0; k < 3; k++) residual[i][k] -= flux[i - 1][k];
            } else {
                // Interface: the oxide is a capacitor between the gate and the interface, the fixed
                // oxide charge sits here, and no carrier flux enters the oxide
                residual[i][0] += Q_f - C_ox * (V(gi, time) - V_gate);
                local[i][0][0] += C_ox;
            }

            local[i][0][0] += q * (dn + dp);
            local[i][0][1] = - q * dn;
            local[i][0][2] = - q * dp;

            local[i][1][0] = - ((r_n + 1 / dt) * dn - r_p * dp);
            local[i][1][1] = (r_n + 1 / dt) * dn + 1e-300;  // 1e-300 keeps the system solvable if n underflows to zero
            local[i][1][2] = - r_p * dp;

            local[i][2][0] = r_n * dn - (r_p + 1 / dt) * dp;
            local[i][2][1] = - r_n * dn;
            local[i][2][2] = (r_p + 1 / dt) * dp + 1e-300;
        }

        std::vector<Vec3<double>> delta = solveBlockFluxSystem(up, down, local, residual);

        double max_change = 0;
        for (int i = 0; i < N - 1; i++){
            int gi = i + start_SC;  // global index
            for (int k = 0; k < 3; k++){
                double change = delta[i][k];
                if (std::isnan(change)){
                    std::cerr << "transientStep: NaN at node " << gi << std::endl;
                    return -1;
                }
                // Damp large steps so that the exponentials in n, p don't overshoot
                if (fabs(change) > 5 * Vt) change = (change > 0 ? 5 : -5) * Vt;
                max_change = std::max(max_change, fabs(change));
                if (k == 0) V(gi, time) += change;
                else Phi(gi, k - 1) += change;
            }
        }
        converged = max_change < tol_newton;
        iter++;
    }

    // Store the new state. In the oxide the potential falls linearly from the gate to the interface.
    for (int i = 0; i < N - 1; i++){
        int gi = i + start_SC;  // global index
        n(gi, time) = ni * exp((V(gi, time) - Phi(gi, 0)) / Vt);
        p(gi, time) = ni * exp((Phi(gi, 1) - V(gi, time)) / Vt);
    }
    for (int i = 0; i < start_SC; i++){
        double fraction = (x[i] - x[0]) / (x[start_SC] - x[0]);
        V(i, time) = V_gate + (V(start_SC, time) - V_gate) * fraction;
    }

    return converged ? iter : -1;
}
