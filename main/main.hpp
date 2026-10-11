#ifndef MAIN_HPP  // Header guard to prevent multiple inclusions
#define MAIN_HPP

#include "../math/math.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>

// Physical constants
constexpr double q = 1.6e-19;     // Elementary charge
constexpr double kb = 1.38e-23;  // Boltzmann constant
constexpr double e_0 = 8.85e-12; // Vacuum permittivity
constexpr double e_si = 11.7 * e_0; // Silicon permittivity
constexpr double e_ox = 3.9 * e_0; // Oxide permittivity
constexpr double chi_si = 4.05; // Electron affinity of silicon (V)
constexpr double E_g = 1.12;    // Band gap of silicon (eV)
constexpr double pi = 3.14159265358979323846;

// Numerical constants of the steady-state solver
constexpr double tol_np = 1e-8; // Relative tolerance on n, p
constexpr double tol_V = 1e-10; // Absolute tolerance on V (V)
constexpr int max_iter = 1000;
constexpr double max_step = 0.1; // Largest change of gate voltage between two steady-state solutions (V)

// Simulation parameters. They are set, and described, in settings.cpp.
extern double NA, ND, phi_m, Q_f, t_ox, T, Vt, ni, mu_n, mu_p, tau_n, tau_p;   // device
extern int x_points, start_SC;                                                 // mesh
extern double a_x, b_x;
extern std::vector<double> x;
extern double VG_l, VG_u;                                                      // bias sweep
extern int n_VG;
extern std::vector<double> freqs;
extern std::vector<double> sweep_VG;                                           // frequency sweep
extern double f_min, f_max;
extern int points_per_decade;
extern std::vector<double> sweep_rates;                                        // transient
extern int substeps;
extern std::vector<double> step_fractions;
extern double t_ramp, t_first, t_last, growth;

// Conditions that can be solved. Each one keeps its solution in its own column of V, n, p:
// column 0 = equilibrium (VG = 0)
// column 1 = DC bias (steady state)
// column 2 = deep depletion: DC bias swept too fast for the minority carrier to follow, so the
//            minority carrier keeps the profile it had at the start of the sweep
// column 3 = transient: the state at the present time of a time-dependent simulation
enum Condition {
    EQUILIBRIUM,
    DC,
    DEEP_DEPLETION,
    TRANSIENT
};

// Column of V, n, p that holds the solution of a condition
inline int column(Condition cond){
    switch (cond) {
        case DC: return 1;
        case DEEP_DEPLETION: return 2;
        case TRANSIENT: return 3;
        default: return 0;
    }
}

// Column that keeps the steady state at the accumulation end of the gate-voltage range,
// from which the deep-depletion and transient sweeps start
constexpr int START_COLUMN = 4;
constexpr int n_columns = 5;

// Width of semiconductor owned by node i (half of each neighbouring interval). Zero inside the oxide.
double controlVolume(const std::vector<double>& x, int i);

// Solves Poisson's equation and the carrier continuity equations of a steady-state condition
// (EQUILIBRIUM, DC or DEEP_DEPLETION) one after the other until they agree
void consistent(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond);

// Moves the gate node to potential V_gate in steps of at most max_step, each one starting from the last solution
void rampGate(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond, double V_gate);

// Output (output.cpp)
// V, n, p of column `time_step` at every grid node
void save_x_based(const std::string& filename, const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, int time_step);
// Real and imaginary parts of the AC phasors at every grid node
void save_ac(const std::string& filename, const std::vector<double>& x, const std::vector<std::complex<double>>& V_ac,
             const std::vector<std::complex<double>>& n_ac, const std::vector<std::complex<double>>& p_ac);
// Table: each row of data is the value of `first_label` followed by one value for each label
void save_table(const std::string& filename, const std::string& first_label, const std::vector<std::string>& labels,
                const std::vector<std::vector<double>>& data);

#endif
