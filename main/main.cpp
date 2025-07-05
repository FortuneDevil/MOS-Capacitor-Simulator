#include "main.hpp"

// Define global variables (defined in main.cpp, declared as extern in main.hpp)
double NA = 1e16, ND = 1e16; // Example doping concentrations (adjust as needed)
double T = 300;  // Temperature in Kelvin (default)
double Vt = kb * T / q;  // Thermal voltage (calculated from constants)
int x_points = 500;  // Number of grid points
int start_SC = 10;   // Starting point for semiconductor region
std::vector<double> x(x_points, 0);  // Initialize x, all set to 0

void save_to_csv(const std::string& filename, const std::vector<double>& x, Matrix& V, Matrix& n, Matrix& p, int time_step) {
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

int main(){
    double freq = 1e3;  // Frequency in rad/s
    double delta_t = 1e-3/freq;
}