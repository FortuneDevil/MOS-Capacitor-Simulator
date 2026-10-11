#ifndef MAIN_HPP  // Header guard to prevent multiple inclusions
#define MAIN_HPP

#include "../math/math.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <complex>

// Constants
constexpr double a_x = 1e-10;   // First mesh spacing in the semiconductor (m)
constexpr double b_x = 1.01;    // Ratio between consecutive mesh spacings
constexpr double tol_np = 1e-8; // Relative tolerance on n, p
constexpr double tol_V = 1e-10; // Absolute tolerance on V (V)
constexpr int max_iter = 1000;
constexpr double q = 1.6e-19;     // Elementary charge
constexpr double kb = 1.38e-23;  // Boltzmann constant
constexpr double e_0 = 8.85e-12; // Vacuum permittivity
constexpr double e_si = 11.7 * e_0; // Silicon permittivity
constexpr double e_ox = 3.9 * e_0; // Oxide permittivity
constexpr double tau_n = 1e-6;  // Electron lifetime
constexpr double tau_p = 1e-6;  // Hole lifetime
constexpr double pi = 3.14159265358979323846;
constexpr double chi_si = 4.05; // Electron affinity of silicon (V)
constexpr double E_g = 1.12;    // Band gap of silicon (eV)
constexpr double max_step = 0.1; // Largest change of gate voltage between two steady-state solutions (V)

// Global variables
extern double NA, ND;          // Doping: NA > ND gives a p-type substrate, ND > NA an n-type one
extern double phi_m;           // Gate work function (V)
extern double Q_f;             // Fixed oxide charge at the interface (C/m^2)
extern double T;               // Temperature in Kelvin (default 300)
extern double Vt;              // Thermal voltage
extern int x_points;           // Number of grid intervals (grid nodes are 0..x_points)
extern int n_VG;               // Number of points for VG
extern int start_SC;           // Node at the oxide-semiconductor interface
extern std::vector<double> x;  // Grid of positions
extern double ni;
extern double mu_n, mu_p;

// Enum for the steady-state conditions. They select the column of V, n, p that is solved:
// column 0 = equilibrium (VG = 0), column 1 = DC bias
enum Condition {
    EQUILIBRIUM,
    DC
};

// Width of semiconductor owned by node i (half of each neighbouring interval). Zero inside the oxide.
double controlVolume(const std::vector<double>& x, int i);

void consistent(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond);

// Moves the gate node to potential V_gate in steps of at most max_step, each one starting from the last solution
void rampGate(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond, double V_gate);

void save_x_based(const std::string& filename, const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, int time_step);
void save_ac(const std::string& filename, const std::vector<double>& x, const std::vector<std::complex<double>>& V_ac,
             const std::vector<std::complex<double>>& n_ac, const std::vector<std::complex<double>>& p_ac);
// Bias sweep: one row per VG, one column of `quantity` per frequency
void save_capacitance_data(const std::string& filename, const std::string& quantity, const std::vector<double>& freqs,
                           const std::vector<std::vector<double>>& data);
// Frequency sweep: one row per frequency, columns C and G/w for each VG
void save_frequency_sweep(const std::string& filename, const std::vector<double>& biases, const std::vector<std::vector<double>>& data);

#endif
