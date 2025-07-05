#ifndef MAIN_HPP  // Header guard to prevent multiple inclusions
#define MAIN_HPP

#include "../math/math.hpp"
#include <fstream>
#include <vector>
#include <string>

// Constants
constexpr double a_x = 1e-10;
constexpr double b_x = 1.01;
constexpr double tol_np = 1e-3;
constexpr double tol_V = 1e-6;
constexpr int max_iter = 10000;
constexpr double q = 1.6e-19;     // Elementary charge
constexpr double kb = 1.38e-23;  // Boltzmann constant
constexpr double e_0 = 8.85e-14; // Vacuum permittivity
constexpr double e_si = 11.7 * e_0; // Silicon permittivity
constexpr double e_ox = 3.9 * e_0; // Oxide permittivity
constexpr double tau_n = 1e-6;  // Electron lifetime
constexpr double tau_p = 1e-6;  // Hole lifetime

// Global variables
extern double NA, ND;
extern double T;               // Temperature in Kelvin (default 300)
extern double Vt;              // Thermal voltage
extern int x_points;           // Number of grid points
extern int start_SC;           // Starting point for semiconductor region
extern std::vector<double> x;  // Grid of positions

// Enum for simulation conditions
enum Condition {
    EQUILIBRIUM,
    DC,
    AC
};

void save_to_csv(const std::string& filename, const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, int time_step);

#endif