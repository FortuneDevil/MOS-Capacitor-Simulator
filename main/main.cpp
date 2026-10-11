#include "main.hpp"
#include "../poisson/poisson.hpp"
#include "../carrier/carrier.hpp"

// Define global variables (defined in main.cpp, declared as extern in main.hpp)
double NA = 1e22, ND = 1e16; // Example doping concentrations (adjust as needed)
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

void save_capacitance_data(const std::string& filename, const std::vector<double>& freqs, const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    // One capacitance column per frequency
    file << "VG";
    for (double f : freqs) file << ",C (f = " << f << " Hz)";
    file << "\n";
    for (const auto& row : data) {
        file << std::scientific << std::setprecision(6) << row[0];
        for (size_t k = 1; k < row.size(); k++) file << "," << row[k];
        file << "\n";
    }

    file.close();
    std::cout << "Capacitance data saved to " << filename << std::endl;
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
    std::cout << "EQUILIBRIUM\n\n";
    consistent(x, V, n, p, EQUILIBRIUM); // solve for equilibrium
    save_x_based("csv_files/eq.csv", x, V, n, p, 0);

    // Start the DC solution from equilibrium
    V.copyColumn(0, 1); // Copies column 0 to 1
    n.copyColumn(0, 1); // Copies column 0 to 1
    p.copyColumn(0, 1); // Copies column 0 to 1

    std::vector<std::vector<double>> vg_c_data;
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    // Now apply a DC bias
    double VG, VG_u = 4, VG_l = -2; // VG ranges from VG_l to VG_u
    double max_step = 0.1;  // largest change of gate voltage between two DC solutions
    for (int i = 0; i < n_VG; i++){
        VG = VG_l + (VG_u - VG_l) * i / (n_VG - 1);
        std::cout << "V_bias = " << VG << std::endl;

        // Move the gate from the previous bias to VG in small steps, each one starting from the last solution
        int steps = std::max(1, (int)ceil(fabs(VG - V(0, 1)) / max_step));
        double step = (VG - V(0, 1)) / steps;
        for (int k = 1; k <= steps; k++){
            V(0, 1) = (k == steps) ? VG : V(0, 1) + step;
            consistent(x, V, n, p, DC);
        }

        // Small-signal AC analysis around this DC solution
        std::vector<double> row = {VG};
        for (double f : freqs){
            std::complex<double> C = smallSignalAC(x, V, n, p, 2 * pi * f, V_ac, n_ac, p_ac);
            row.push_back(C.real());
        }
        vg_c_data.push_back(row);
    }
    save_x_based("csv_files/dc.csv", x, V, n, p, 1);
    save_ac("csv_files/ac.csv", x, V_ac, n_ac, p_ac);
    save_capacitance_data("csv_files/VG_C.csv", freqs, vg_c_data);
    return 0;
}
