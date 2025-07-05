#include "main.hpp"
#include "../poisson/poisson.hpp"
#include "../carrier/carrier.hpp"

// Define global variables (defined in main.cpp, declared as extern in main.hpp)
double NA = 1e16, ND = 1e16; // Example doping concentrations (adjust as needed)
double T = 300;  // Temperature in Kelvin (default)
double Vt = kb * T / q;  // Thermal voltage (calculated from constants)
int x_points = 500;  // Number of grid points
int start_SC = 10;   // Starting point for semiconductor region
std::vector<double> x(x_points, 0);  // Initialize x, all set to 0

void consistent(const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, Condition cond){
    int iter = 0;
    while (iter < max_iter){
        Matrix temp_V = V, temp_n = n, temp_p = p;
        poissonSolver(x, V, n, p, cond);
        carrier(x, V, n, p, cond);
        if (max_relative_error(V, temp_V) < tol_V &&
            max_relative_error(n, temp_n) < tol_np &&
            max_relative_error(p, temp_p) < tol_np){
                std::cout << "Converged at iteration " << iter << std::endl;
                break;
        }
        iter++;
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

void save_capacitance_data(const std::string& filename, const std::vector<std::pair<double, double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    file << "VG,C\n";
    for (const auto& [vg, c] : data) {
        file << std::scientific << std::setprecision(6)
             << vg << "," << c << "\n";
    }

    file.close();
    std::cout << "Capacitance data saved to " << filename << std::endl;
}

double ni = 1.5e10;

int main(){
    double freq = 1e3;  // Frequency in rad/s
    double delta_t = 1e-3/freq;
    double t_ox = 10e-9;    //10nm oxide
    Matrix V(x_points + 1, 3), n(x_points + 1, 3), p(x_points + 1, 3);
    for (int i = 0; i <= start_SC; i++){
        x[i] = i * t_ox / start_SC;
    }
    for (int i = start_SC + 1; i <= x_points; i++){
        x[i] = t_ox + a_x * (pow(b_x, (i - start_SC) - 1)) / (b_x - 1);
        n(i, 0) = pow(ni, 2) / NA;
        p(i, 0) = NA;
    }
    consistent(x, V, n, p, EQUILIBRIUM); // solve for equilibrium
    
    std::vector<std::pair<double, double>> vg_c_data;
    // Now apply a DC bias
    double VG, VG_u = 2, VG_l = -2; // VG ranges from VG_l to VG_u
    for (int i = 0; i < n_VG; i++){
        VG = VG_l + (VG_u - VG_l) * i / (n_VG - 1);
        V.copyColumn(0, 1); // Copies column 0 to 1
        n.copyColumn(0, 1); // Copies column 0 to 1
        p.copyColumn(0, 1); // Copies column 0 to 1
        
        V(0, 1) = VG;
        consistent(x, V, n, p, DC);

        V.copyColumn(1, 2); // Copies column 1 to 2
        n.copyColumn(1, 2); // Copies column 1 to 2
        p.copyColumn(1, 2); // Copies column 1 to 2
        
        V(0, 2) += V(0, 2)/1e3;    // add a small ac
        consistent(x, V, n, p, AC);

        // Calculate the capacitance
        double C = e_ox / (x[1] - x[0]) * (1 - (V(1, 2) - V(1, 1)) / (V(0, 2) - V(0, 1)));

        vg_c_data.emplace_back(VG, C);
    }
    save_x_based("../csv_files/eq.csv", x, V, n, p, 0);
    save_x_based("../csv_files/dc.csv", x, V, n, p, 1);
    save_x_based("../csv_files/ac.csv", x, V, n, p, 2);
    save_capacitance_data("../csv_files/VG_C.csv", vg_c_data);
    return 0;
}