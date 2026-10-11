#include "main.hpp"
#include "../poisson/poisson.hpp"
#include "../carrier/carrier.hpp"

// Define global variables (defined in main.cpp, declared as extern in main.hpp)
// Doping concentrations (adjust as needed). NA > ND gives a p-type substrate, ND > NA an n-type one.
double NA = 1e22, ND = 1e16;
double phi_m = 4.1;         // Gate work function (V). 4.1 is aluminium.
double Q_f = q * 1e15;      // Fixed oxide charge at the interface (C/m^2). Here 1e15 charges/m^2.
double T = 300;  // Temperature in Kelvin (default)
double Vt = kb * T / q;  // Thermal voltage (calculated from constants)
int x_points = 500;  // Number of grid intervals
int start_SC = 10;   // Node at the oxide-semiconductor interface
std::vector<double> x(x_points + 1, 0);  // Initialize x, all set to 0 (indices 0..x_points)

int n_VG = 81;

double ni = 1.5e10*1e6;
double mu_n = 0.1360;
double mu_p = 0.0480;

double controlVolume(const std::vector<double>& x, int i){
    int last = x.size() - 1;
    if (i < start_SC) return 0;
    if (i == start_SC) return (x[i + 1] - x[i]) / 2;
    if (i == last) return (x[i] - x[i - 1]) / 2;
    return (x[i + 1] - x[i - 1]) / 2;
}

void consistent(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond){
    int iter = 0;
    while (iter < max_iter){
        double change_V = poissonSolver(x, V, n, p, cond);
        double change_np = carrier(x, V, n, p, cond);
        if (change_V < tol_V && change_np < tol_np){
            std::cout << "Converged at iteration " << iter << std::endl;
            return;
        }
        iter++;
    }
    std::cerr << "WARNING: not converged after " << max_iter << " iterations" << std::endl;
}

void rampGate(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond, double V_gate){
    int time = (cond == DC) ? 1 : 0;
    int steps = std::max(1, (int)ceil(fabs(V_gate - V(0, time)) / max_step));
    double step = (V_gate - V(0, time)) / steps;
    for (int k = 1; k <= steps; k++){
        V(0, time) = (k == steps) ? V_gate : V(0, time) + step;
        consistent(x, V, n, p, cond);
    }
}

void save_x_based(const std::string& filename, const std::vector<double>& x, const Matrix& V, const Matrix& n, const Matrix& p, int time_step) {
    // Open the CSV file for writing
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Write the header to the CSV file
    file << "x,V,n,p\n";

    // Write data for each x[i], V[i], n[i], p[i]
    int distance = x.size();  // Total number of points in x
    for (int i = 0; i < distance; i++) {
        // Get corresponding values from V, n, and p matrices at the given time_step
        double voltage = (i < V.size()) ? V(i, time_step) : 0.0;
        double electron = (i < n.size()) ? n(i, time_step) : 0.0;
        double hole = (i < p.size()) ? p(i, time_step) : 0.0;

        // Write the row: x, V, n, p values
        file << std::scientific << std::setprecision(6)  // Print in scientific notation
             << x[i] << ","
             << voltage << ","
             << electron << ","
             << hole << "\n";
    }

    // Close the file after writing
    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

void save_ac(const std::string& filename, const std::vector<double>& x, const std::vector<std::complex<double>>& V_ac,
             const std::vector<std::complex<double>>& n_ac, const std::vector<std::complex<double>>& p_ac) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    // Real and imaginary parts of the AC phasors, per volt of gate signal
    file << "x,V_re,V_im,n_re,n_im,p_re,p_im\n";
    for (size_t i = 0; i < x.size(); i++) {
        file << std::scientific << std::setprecision(6)
             << x[i] << ","
             << V_ac[i].real() << "," << V_ac[i].imag() << ","
             << n_ac[i].real() << "," << n_ac[i].imag() << ","
             << p_ac[i].real() << "," << p_ac[i].imag() << "\n";
    }

    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}

void save_capacitance_data(const std::string& filename, const std::string& quantity, const std::vector<double>& freqs,
                           const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    // One column per frequency
    file << "VG";
    for (double f : freqs) file << "," << quantity << " (f = " << f << " Hz)";
    file << "\n";
    for (const auto& row : data) {
        file << std::scientific << std::setprecision(6) << row[0];
        for (size_t k = 1; k < row.size(); k++) file << "," << row[k];
        file << "\n";
    }

    file.close();
    std::cout << quantity << " data saved to " << filename << std::endl;
}

void save_frequency_sweep(const std::string& filename, const std::vector<double>& biases, const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    // Two columns per bias: capacitance and conductance / omega
    file << "f";
    for (double VG : biases) file << ",C (VG = " << VG << " V),G/w (VG = " << VG << " V)";
    file << "\n";
    for (const auto& row : data) {
        file << std::scientific << std::setprecision(6) << row[0];
        for (size_t k = 1; k < row.size(); k++) file << "," << row[k];
        file << "\n";
    }

    file.close();
    std::cout << "Frequency sweep saved to " << filename << std::endl;
}

int main(){
    // Frequencies (Hz) of the small-signal analysis. The last one is the one saved in ac.csv.
    std::vector<double> freqs = {1e-3, 1e-2, 5e-2, 1e-1, 1, 10, 1e3, 1e6};
    double t_ox = 10e-9;    //10nm oxide
    Matrix V(x_points + 1, 2), n(x_points + 1, 2), p(x_points + 1, 2);

    // Charge-neutral bulk: p - n = NA - ND and n p = ni^2
    double doping = NA - ND;
    double p_bulk, n_bulk;
    if (doping >= 0){
        p_bulk = (doping + sqrt(doping * doping + 4 * ni * ni)) / 2;
        n_bulk = ni * ni / p_bulk;
    } else {
        n_bulk = (-doping + sqrt(doping * doping + 4 * ni * ni)) / 2;
        p_bulk = ni * ni / n_bulk;
    }

    // Mesh: uniform in the oxide, spacing a_x * b_x^k in the semiconductor
    for (int i = 0; i <= start_SC; i++){
        x[i] = i * t_ox / start_SC;
    }
    for (int i = start_SC; i <= x_points; i++){
        if (i > start_SC) x[i] = t_ox + a_x * (pow(b_x, i - start_SC) - 1) / (b_x - 1);
        n(i, 0) = n_bulk;
        p(i, 0) = p_bulk;
    }

    // Work function difference between gate and semiconductor. The semiconductor work function is
    // chi + (Ec - Ef)/q, with the Fermi level set by the bulk doping.
    double phi_s = chi_si + E_g / 2 - Vt * log(n_bulk / ni);
    double phi_ms = phi_m - phi_s;
    // V is measured from the bulk, so a gate voltage VG puts the gate node at VG - phi_ms
    std::cout << "phi_ms = " << phi_ms << " V, flat-band voltage = " << phi_ms - Q_f * t_ox / e_ox << " V\n";

    std::cout << "EQUILIBRIUM\n\n";
    rampGate(x, V, n, p, EQUILIBRIUM, 0 - phi_ms); // solve for equilibrium (VG = 0)
    save_x_based("csv_files/eq.csv", x, V, n, p, 0);

    // Start the DC solution from equilibrium
    V.copyColumn(0, 1); // Copies column 0 to 1
    n.copyColumn(0, 1); // Copies column 0 to 1
    p.copyColumn(0, 1); // Copies column 0 to 1

    // Capacitance C and conductance over omega G/w versus VG, one column per frequency
    std::vector<std::vector<double>> vg_c_data, vg_g_data;
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    // Now apply a DC bias
    double VG, VG_u = 4, VG_l = -4; // VG ranges from VG_l to VG_u
    for (int i = 0; i < n_VG; i++){
        VG = VG_l + (VG_u - VG_l) * i / (n_VG - 1);
        std::cout << "V_bias = " << VG << std::endl;
        rampGate(x, V, n, p, DC, VG - phi_ms);

        // Small-signal AC analysis around this DC solution.
        // The admittance is j w C_complex, so C = Re(C_complex) and G/w = -Im(C_complex).
        std::vector<double> row_c = {VG}, row_g = {VG};
        for (double f : freqs){
            std::complex<double> C = smallSignalAC(x, V, n, p, 2 * pi * f, V_ac, n_ac, p_ac);
            row_c.push_back(C.real());
            row_g.push_back(-C.imag());
        }
        vg_c_data.push_back(row_c);
        vg_g_data.push_back(row_g);
    }
    save_x_based("csv_files/dc.csv", x, V, n, p, 1);
    save_ac("csv_files/ac.csv", x, V_ac, n_ac, p_ac);
    save_capacitance_data("csv_files/VG_C.csv", "C", freqs, vg_c_data);
    save_capacitance_data("csv_files/VG_G.csv", "G/w", freqs, vg_g_data);

    // Frequency sweep at a few fixed gate voltages
    std::vector<double> sweep_VG = {-2, -0.5, 0, 4};
    double f_min = 1e-4, f_max = 1e7;   // Hz
    int points_per_decade = 10;
    int n_f = (int)round(log10(f_max / f_min) * points_per_decade) + 1;
    std::vector<std::vector<double>> f_data(n_f);
    for (int k = 0; k < n_f; k++){
        f_data[k] = {f_min * pow(10.0, (double)k / points_per_decade)};
    }
    for (double bias : sweep_VG){
        std::cout << "Frequency sweep at V_bias = " << bias << std::endl;
        rampGate(x, V, n, p, DC, bias - phi_ms);
        for (int k = 0; k < n_f; k++){
            std::complex<double> C = smallSignalAC(x, V, n, p, 2 * pi * f_data[k][0], V_ac, n_ac, p_ac);
            f_data[k].push_back(C.real());
            f_data[k].push_back(-C.imag());
        }
    }
    save_frequency_sweep("csv_files/f_C.csv", sweep_VG, f_data);
    return 0;
}
