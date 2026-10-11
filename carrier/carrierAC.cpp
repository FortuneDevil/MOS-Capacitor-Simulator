#include "carrier.hpp"
#include <array>

// Extended precision: the AC fluxes are many orders of magnitude larger than the j*omega terms they balance
typedef long double Real;
typedef std::complex<Real> Complex;
typedef std::array<Complex, 3> Vec3;                // (V_ac, n_ac, p_ac) of one node
typedef std::array<std::array<Complex, 3>, 3> Mat3; // rows: Poisson, electron continuity, hole continuity

static Mat3 zeroMat(){
    Mat3 M;
    for (auto& row : M) row.fill(Complex(0));
    return M;
}

static Mat3 operator+(const Mat3& A, const Mat3& B){
    Mat3 C;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) C[i][j] = A[i][j] + B[i][j];
    return C;
}

static Mat3 operator*(const Mat3& A, const Mat3& B){
    Mat3 C = zeroMat();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
}

static Vec3 operator*(const Mat3& A, const Vec3& v){
    Vec3 w;
    for (int i = 0; i < 3; i++) w[i] = A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2];
    return w;
}

// Inverse by Gauss-Jordan elimination with scaled partial pivoting
static Mat3 inverse(Mat3 A){
    Mat3 inv = zeroMat();
    Real scale[3];
    for (int i = 0; i < 3; i++){
        inv[i][i] = 1;
        scale[i] = std::max({std::abs(A[i][0]), std::abs(A[i][1]), std::abs(A[i][2])});
    }
    for (int i = 0; i < 3; i++){
        int max_row = i;
        for (int r = i + 1; r < 3; r++){
            if (std::abs(A[r][i]) / scale[r] > std::abs(A[max_row][i]) / scale[max_row]) max_row = r;
        }
        std::swap(A[i], A[max_row]);
        std::swap(inv[i], inv[max_row]);
        std::swap(scale[i], scale[max_row]);

        Complex pivot = A[i][i];
        for (int j = 0; j < 3; j++){
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for (int r = 0; r < 3; r++){
            if (r == i) continue;
            Complex factor = A[r][i];
            for (int j = 0; j < 3; j++){
                A[r][j] -= factor * A[i][j];
                inv[r][j] -= factor * inv[i][j];
            }
        }
    }
    return inv;
}

std::complex<double> smallSignalAC(const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, double omega,
                                   std::vector<std::complex<double>>& V_ac, std::vector<std::complex<double>>& n_ac,
                                   std::vector<std::complex<double>>& p_ac){
    int N = V.size() - start_SC;    // semiconductor nodes, from the interface to the bulk contact
    const Real vt = Vt;
    const Complex jw(0, omega);

    // With X_i = (V_ac, n_ac, p_ac) at node i, the three linearised equations of node i read
    //   F(i) - F(i-1) - local[i] * X_i = 0,   F(i) = up[i] * X_(i+1) - down[i] * X_i
    // where F(i) holds the AC fluxes (displacement, electron, hole) between nodes i and i + 1.
    std::vector<Mat3> up(N, zeroMat()), down(N, zeroMat()), local(N, zeroMat());

    for (int i = 0; i < N - 1; i++){
        int gi = i + start_SC;  // global index
        Real delta_x_i = x[gi + 1] - x[gi];
        Real dphi = ((Real)V(gi + 1, 1) - (Real)V(gi, 1)) / vt;
        Real B_plus = bernoulli(dphi), B_minus = bernoulli(-dphi);
        Real dB_plus = dBernoulli(dphi), dB_minus = dBernoulli(-dphi);

        // Poisson: e_si * (V_(i+1) - V_i) / delta_x
        up[i][0][0] = down[i][0][0] = (Real)e_si / delta_x_i;

        // Electron flux: nfactor * [B(dphi) n_(i+1) - B(-dphi) n_i + g_n (V_(i+1) - V_i) / Vt]
        Real nfactor = vt * mu_n / delta_x_i;
        Real g_n = dB_plus * n(gi + 1, 1) + dB_minus * n(gi, 1);
        up[i][1][1] = nfactor * B_plus;
        down[i][1][1] = nfactor * B_minus;
        up[i][1][0] = down[i][1][0] = nfactor * g_n / vt;

        // Hole flux (taken towards node i): pfactor * [B(-dphi) p_(i+1) - B(dphi) p_i - g_p (V_(i+1) - V_i) / Vt]
        Real pfactor = vt * mu_p / delta_x_i;
        Real g_p = dB_plus * p(gi, 1) + dB_minus * p(gi + 1, 1);
        up[i][2][2] = pfactor * B_minus;
        down[i][2][2] = pfactor * B_plus;
        up[i][2][0] = down[i][2][0] = - pfactor * g_p / vt;

        // Charge in Poisson's equation
        Real h = controlVolume(x, gi);
        local[i][0][1] = q * h;
        local[i][0][2] = - q * h;

        // Linearised SRH recombination R_ac = r_n n_ac + r_p p_ac, and the j*omega storage term
        Real n0 = n(gi, 1), p0 = p(gi, 1);
        Real srh = tau_p * (n0 + ni) + tau_n * (p0 + ni);
        Real R0 = (n0 * p0 - (Real)ni * ni) / srh;
        Real r_n = (p0 - R0 * tau_p) / srh;
        Real r_p = (n0 - R0 * tau_n) / srh;
        local[i][1][1] = h * (r_n + jw);
        local[i][1][2] = h * r_p;
        local[i][2][1] = h * r_n;
        local[i][2][2] = h * (r_p + jw);
    }

    // Forward Elimination, never subtracting the large fluxes from each other.
    // After it, node i obeys: up[i] * X_(i+1) - (down[i] + E_i) * X_i = -beta_i
    std::vector<Mat3> pivot_inv(N, zeroMat());
    std::vector<Vec3> beta(N, Vec3{Complex(0), Complex(0), Complex(0)});

    // Interface node: the oxide is a capacitor between the gate (1 V of AC signal) and the interface.
    // No carrier flux enters the oxide.
    Real C_ox = e_ox / (x[start_SC] - x[0]);
    Mat3 E = local[0];
    E[0][0] += C_ox;
    beta[0][0] = C_ox;

    for (int i = 0; i < N - 1; i++){
        pivot_inv[i] = inverse(down[i] + E);
        if (i + 1 < N - 1){
            E = local[i + 1] + E * pivot_inv[i] * up[i];
            beta[i + 1] = down[i] * (pivot_inv[i] * beta[i]);
        }
    }

    // Backward Substitution. The bulk contact is AC ground with no change in n, p.
    std::vector<Vec3> X(N, Vec3{Complex(0), Complex(0), Complex(0)});
    for (int i = N - 2; i >= 0; i--){
        Vec3 right = up[i] * X[i + 1];
        for (int k = 0; k < 3; k++) right[k] += beta[i][k];
        X[i] = pivot_inv[i] * right;
    }

    auto toDouble = [](const Complex& z){ return std::complex<double>(z.real(), z.imag()); };
    V_ac.assign(V.size(), 0); n_ac.assign(V.size(), 0); p_ac.assign(V.size(), 0);
    for (int i = 0; i < N; i++){
        int gi = i + start_SC;  // global index
        V_ac[gi] = toDouble(X[i][0]);
        n_ac[gi] = toDouble(X[i][1]);
        p_ac[gi] = toDouble(X[i][2]);
    }
    // In the oxide the AC potential falls linearly from the gate to the interface
    for (int i = 0; i < start_SC; i++){
        Real fraction = (x[i] - x[0]) / (x[start_SC] - x[0]);
        V_ac[i] = toDouble(Complex(1) + (X[0][0] - Complex(1)) * fraction);
    }

    // Gate charge is the charge on the oxide capacitor, and the gate signal is 1 V
    return toDouble(C_ox * (Complex(1) - X[0][0]));
}
