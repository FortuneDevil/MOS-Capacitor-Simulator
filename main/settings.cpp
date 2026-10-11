#include "main.hpp"

// All the parameters of a simulation are set in this file. SI units throughout.

// ---------------- Device ----------------
// Doping concentrations (m^-3). NA > ND gives a p-type substrate, ND > NA an n-type one.
double NA = 1e22, ND = 1e16;
double phi_m = 4.1;         // Gate work function (V). 4.1 is aluminium.
double Q_f = q * 1e15;      // Fixed oxide charge at the interface (C/m^2). Here 1e15 charges/m^2.
double t_ox = 10e-9;        // Oxide thickness (m)
double T = 300;             // Temperature (K)
double Vt = kb * T / q;     // Thermal voltage (V)
double ni = 1.5e10*1e6;     // Intrinsic carrier concentration (m^-3)
double mu_n = 0.1360;       // Electron mobility (m^2/Vs)
double mu_p = 0.0480;       // Hole mobility (m^2/Vs)
double tau_n = 1e-6;        // Electron lifetime (s)
double tau_p = 1e-6;        // Hole lifetime (s)

// ---------------- Mesh ----------------
int x_points = 500;         // Number of grid intervals (grid nodes are 0..x_points)
int start_SC = 10;          // Node at the oxide-semiconductor interface (the oxide has start_SC intervals)
double a_x = 1e-10;         // First mesh spacing in the semiconductor (m)
double b_x = 1.01;          // Ratio between consecutive mesh spacings
std::vector<double> x(x_points + 1, 0);  // Grid of positions, filled in main.cpp

// ---------------- Bias sweep (C-V, G-V, deep depletion) ----------------
double VG_l = -4, VG_u = 4; // Gate voltage ranges from VG_l to VG_u (V)
int n_VG = 81;              // Number of gate voltages
// Frequencies (Hz) of the small-signal analysis. The last one is the "high frequency": it is the one
// saved in ac.csv and the one used for the deep-depletion and transient capacitances.
std::vector<double> freqs = {1e-3, 1e-2, 5e-2, 1e-1, 1, 10, 1e3, 1e6};

// ---------------- Frequency sweep ----------------
std::vector<double> sweep_VG = {-2, -0.5, 0, 4};    // Gate voltages of the frequency sweep (V)
double f_min = 1e-4, f_max = 1e7;                   // Frequency range (Hz)
int points_per_decade = 10;

// ---------------- Transient: C-V at a finite sweep rate ----------------
std::vector<double> sweep_rates = {1e-3, 1e-2, 1e-1, 1, 10};  // Gate sweep rates (V/s)
int substeps = 4;           // Time steps between two gate voltages of the sweep

// ---------------- Transient: C-t after a gate step ----------------
// The gate is stepped from the accumulation end of the range to these fractions of the inversion end
std::vector<double> step_fractions = {0.5, 1.0};
double t_ramp = 1e-9;       // Duration of each 0.1 V of the step (s)
double t_first = 1e-6, t_last = 1e3;    // First and last time after the step (s)
double growth = 1.15;       // Ratio between consecutive times
