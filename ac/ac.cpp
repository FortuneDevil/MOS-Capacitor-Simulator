#include "ac.hpp"

// Extended precision, as a safety margin: the AC fluxes are many orders of magnitude larger than the
// j*omega terms they balance (double already gives the same results for the default device)
typedef long double Real;
typedef std::complex<Real> Complex;

std::complex<double> smallSignalAC(const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, int time, double omega,
                                   std::vector<std::complex<double>>& V_ac, std::vector<std::complex<double>>& n_ac,
                                   std::vector<std::complex<double>>& p_ac){
    int N = V.size() - start_SC;    // semiconductor nodes, from the interface to the bulk contact
    const Real vt = Vt;
    const Complex jw(0, omega);

    // The unknowns at node i are X_i = (V_ac, Wn, Wp), where Wn and Wp are the AC quasi-Fermi potentials:
    //   n_ac = n0 (V_ac - Wn) / Vt,   p_ac = p0 (Wp - V_ac) / Vt
    // No DC current flows in a MOS capacitor, so in these variables the drift and diffusion parts of the
    // AC currents combine into a plain conductance times the difference of Wn (or Wp) across each element.
    // The three linearised equations of node i then read
    //   F(i) - F(i-1) - local[i] * X_i = 0,   F(i) = up[i] * X_(i+1) - down[i] * X_i
    // where F(i) holds the AC fluxes (displacement, electron, hole) between nodes i and i + 1.
    std::vector<Mat3<Complex>> up(N, zeroMat3<Complex>()), down(N, zeroMat3<Complex>()), local(N, zeroMat3<Complex>());

    for (int i = 0; i < N - 1; i++){
        int gi = i + start_SC;  // global index
        Real delta_x_i = x[gi + 1] - x[gi];
        Real dphi = ((Real)V(gi + 1, time) - (Real)V(gi, time)) / vt;
        Real B_plus = bernoulli(dphi), B_minus = bernoulli(-dphi);

        // Poisson: e_si * (V_(i+1) - V_i) / delta_x
        up[i][0][0] = down[i][0][0] = (Real)e_si / delta_x_i;

        // Electrons: conductance * (Wn_(i+1) - Wn_i). With zero DC current B(dphi) n_(i+1) = B(-dphi) n_i.
        Real n_face = (B_plus * n(gi + 1, time) + B_minus * n(gi, time)) / 2;
        up[i][1][1] = down[i][1][1] = mu_n * n_face / delta_x_i;

        // Holes: conductance * (Wp_(i+1) - Wp_i). With zero DC current B(-dphi) p_(i+1) = B(dphi) p_i.
        Real p_face = (B_minus * p(gi + 1, time) + B_plus * p(gi, time)) / 2;
        up[i][2][2] = down[i][2][2] = mu_p * p_face / delta_x_i;

        Real h = controlVolume(x, gi);
        Real n0 = n(gi, time), p0 = p(gi, time);
        Real dn = h * n0 / vt, dp = h * p0 / vt;    // h * n_ac = dn * (V_ac - Wn),  h * p_ac = dp * (Wp - V_ac)

        // Charge in Poisson's equation: q h (p_ac - n_ac)
        local[i][0][0] = q * (dn + dp);
        local[i][0][1] = - q * dn;
        local[i][0][2] = - q * dp;

        // Linearised SRH recombination R_ac = r_n n_ac + r_p p_ac, and the j*omega storage term
        Real srh = tau_p * (n0 + ni) + tau_n * (p0 + ni);
        Real R0 = (n0 * p0 - (Real)ni * ni) / srh;
        Real r_n = (p0 - R0 * tau_p) / srh;
        Real r_p = (n0 - R0 * tau_n) / srh;
        // Electrons: + h (R_ac + jw n_ac)
        local[i][1][0] = - ((r_n + jw) * dn - r_p * dp);
        local[i][1][1] = (r_n + jw) * dn;
        local[i][1][2] = - r_p * dp;
        // Holes: - h (R_ac + jw p_ac)
        local[i][2][0] = r_n * dn - (r_p + jw) * dp;
        local[i][2][1] = - r_n * dn;
        local[i][2][2] = (r_p + jw) * dp;

        // Keeps the equations solvable where a carrier concentration underflows to zero
        local[i][1][1] += (Real)1e-300;
        local[i][2][2] += (Real)1e-300;
    }

    // Interface node: the oxide is a capacitor between the gate (1 V of AC signal) and the interface.
    // No carrier flux enters the oxide.
    Real C_ox = e_ox / (x[start_SC] - x[0]);
    std::vector<Vec3<Complex>> source(N, zeroVec3<Complex>());
    local[0][0][0] += C_ox;
    source[0][0] = C_ox;

    // The bulk contact is AC ground with no change in n, p.
    std::vector<Vec3<Complex>> X = solveBlockFluxSystem(up, down, local, source);

    auto toDouble = [](const Complex& z){ return std::complex<double>(z.real(), z.imag()); };
    V_ac.assign(V.size(), 0); n_ac.assign(V.size(), 0); p_ac.assign(V.size(), 0);
    for (int i = 0; i < N; i++){
        int gi = i + start_SC;  // global index
        V_ac[gi] = toDouble(X[i][0]);
        n_ac[gi] = toDouble((Real)n(gi, time) * (X[i][0] - X[i][1]) / vt);
        p_ac[gi] = toDouble((Real)p(gi, time) * (X[i][2] - X[i][0]) / vt);
    }
    // In the oxide the AC potential falls linearly from the gate to the interface
    for (int i = 0; i < start_SC; i++){
        Real fraction = (x[i] - x[0]) / (x[start_SC] - x[0]);
        V_ac[i] = toDouble(Complex(1) + (X[0][0] - Complex(1)) * fraction);
    }

    // Gate charge is the charge on the oxide capacitor, and the gate signal is 1 V
    return toDouble(C_ox * (Complex(1) - X[0][0]));
}
