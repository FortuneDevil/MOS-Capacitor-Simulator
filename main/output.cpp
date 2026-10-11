#include "main.hpp"

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

void save_table(const std::string& filename, const std::string& first_label, const std::vector<std::string>& labels,
                const std::vector<std::vector<double>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing.\n";
        return;
    }

    file << first_label;
    for (const std::string& label : labels) file << "," << label;
    file << "\n";
    for (const auto& row : data) {
        file << std::scientific << std::setprecision(6) << row[0];
        for (size_t k = 1; k < row.size(); k++) file << "," << row[k];
        file << "\n";
    }

    file.close();
    std::cout << "Data saved to " << filename << std::endl;
}
