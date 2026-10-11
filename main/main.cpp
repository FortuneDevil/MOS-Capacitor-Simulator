#include "main.hpp"
#include "../poisson/poisson.hpp"
#include "../carrier/carrier.hpp"
#include "../ac/ac.hpp"
#include "../transient/transient.hpp"

// Quantities derived from the settings
static bool p_type;     // Substrate type. Accumulation is at the lower end of the VG range for p-type, at the upper end for n-type.
static double phi_ms;   // Work function difference. V is measured from the bulk, so a gate voltage VG puts the gate node at VG - phi_ms.

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
    int time = column(cond);
    int steps = std::max(1, (int)ceil(fabs(V_gate - V(0, time)) / max_step));
    double step = (V_gate - V(0, time)) / steps;
    for (int k = 1; k <= steps; k++){
        V(0, time) = (k == steps) ? V_gate : V(0, time) + step;
        consistent(x, V, n, p, cond);
    }
}

// Gate voltage of point i of the bias sweep
static double biasPoint(int i){
    return VG_l + (VG_u - VG_l) * i / (n_VG - 1);
}

// Point of the bias sweep that comes k-th when sweeping from accumulation towards inversion
static int fromAccumulation(int k){
    return p_type ? k : n_VG - 1 - k;
}

// Column label such as "C (f = 1000 Hz)"
static std::string label(const std::string& before, double value, const std::string& after){
    std::ostringstream text;
    text << before << value << after;
    return text.str();
}

static void copyState(Matrix& V, Matrix& n, Matrix& p, int from, int to){
    V.copyColumn(from, to);
    n.copyColumn(from, to);
    p.copyColumn(from, to);
}

// Builds the mesh, fills column 0 with the charge-neutral bulk concentrations and finds phi_ms
static void setupDevice(Matrix& n, Matrix& p){
    p_type = NA > ND;

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
    phi_ms = phi_m - phi_s;
    std::cout << "phi_ms = " << phi_ms << " V, flat-band voltage = " << phi_ms - Q_f * t_ox / e_ox << " V\n";
}

// Steady state at every gate voltage, and small-signal analysis around it at every frequency.
// Fills the tables of capacitance C and conductance over omega G/w (one row per VG, one column per frequency).
static void biasSweep(Matrix& V, Matrix& n, Matrix& p, std::vector<std::vector<double>>& vg_c_data,
                      std::vector<std::vector<double>>& vg_g_data){
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    for (int i = 0; i < n_VG; i++){
        double VG = biasPoint(i);
        std::cout << "V_bias = " << VG << std::endl;
        rampGate(x, V, n, p, DC, VG - phi_ms);

        std::vector<double> row_c = {VG}, row_g = {VG};
        for (double f : freqs){
            std::complex<double> C = smallSignalAC(x, V, n, p, column(DC), 2 * pi * f, V_ac, n_ac, p_ac);
            row_c.push_back(C.real());
            row_g.push_back(-C.imag());
        }
        vg_c_data.push_back(row_c);
        vg_g_data.push_back(row_g);

        // The deep-depletion and transient sweeps start from the steady state at the accumulation end of the range
        if (i == fromAccumulation(0)) copyState(V, n, p, column(DC), START_COLUMN);
    }
    // Profiles at the last gate voltage, and their AC changes at the last frequency
    save_x_based("csv_files/dc.csv", x, V, n, p, column(DC));
    save_ac("csv_files/ac.csv", x, V_ac, n_ac, p_ac);
}

// Deep depletion: the gate is swept from accumulation towards inversion too fast for minority
// carriers to be generated, so no inversion layer forms. Adds one column to the capacitance table.
static void deepDepletionSweep(Matrix& V, Matrix& n, Matrix& p, std::vector<std::vector<double>>& vg_c_data){
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    copyState(V, n, p, START_COLUMN, column(DEEP_DEPLETION));
    for (int k = 0; k < n_VG; k++){
        int i = fromAccumulation(k);
        std::cout << "V_bias = " << biasPoint(i) << std::endl;
        rampGate(x, V, n, p, DEEP_DEPLETION, biasPoint(i) - phi_ms);
        std::complex<double> C = smallSignalAC(x, V, n, p, column(DEEP_DEPLETION), 2 * pi * freqs.back(), V_ac, n_ac, p_ac);
        vg_c_data[i].push_back(C.real());
    }
}

// C and G/w versus frequency at a few fixed gate voltages
static void frequencySweep(Matrix& V, Matrix& n, Matrix& p){
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    int n_f = (int)round(log10(f_max / f_min) * points_per_decade) + 1;
    std::vector<std::vector<double>> f_data(n_f);
    for (int k = 0; k < n_f; k++){
        f_data[k] = {f_min * pow(10.0, (double)k / points_per_decade)};
    }
    std::vector<std::string> labels;
    for (double bias : sweep_VG){
        std::cout << "Frequency sweep at V_bias = " << bias << std::endl;
        rampGate(x, V, n, p, DC, bias - phi_ms);
        for (int k = 0; k < n_f; k++){
            std::complex<double> C = smallSignalAC(x, V, n, p, column(DC), 2 * pi * f_data[k][0], V_ac, n_ac, p_ac);
            f_data[k].push_back(C.real());
            f_data[k].push_back(-C.imag());
        }
        labels.push_back(label("C (VG = ", bias, " V)"));
        labels.push_back(label("G/w (VG = ", bias, " V)"));
    }
    save_table("csv_files/f_C.csv", "f", labels, f_data);
}

// Puts the starting state (steady state at the accumulation end) in the TRANSIENT column
static void startTransient(Matrix& V, Matrix& n, Matrix& p, Matrix& Phi){
    copyState(V, n, p, START_COLUMN, column(TRANSIENT));
    transientInit(V, n, p, Phi);
}

// Moves the gate node linearly to V_gate in `steps` time steps of length dt
static void moveGate(Matrix& V, Matrix& n, Matrix& p, Matrix& Phi, double V_gate, int steps, double dt){
    double from = V(0, column(TRANSIENT));
    for (int s = 1; s <= steps; s++){
        if (transientStep(x, V, n, p, Phi, from + (V_gate - from) * s / steps, dt) < 0){
            std::cerr << "WARNING: transient step not converged" << std::endl;
        }
    }
}

// Transient 1: C-V with the gate swept from accumulation towards inversion at a finite rate.
// The capacitance is the small-signal one at the highest frequency, around the state at each time.
// The equilibrium and deep-depletion curves of the bias sweep are saved alongside for reference.
static void sweepRateCV(Matrix& V, Matrix& n, Matrix& p, Matrix& Phi, const std::vector<std::vector<double>>& vg_c_data){
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    std::vector<std::vector<double>> rate_data(n_VG);
    std::vector<std::string> labels;

    labels.push_back(label("C equilibrium (f = ", freqs.back(), " Hz)"));
    for (int i = 0; i < n_VG; i++){
        rate_data[i] = {vg_c_data[i][0], vg_c_data[i][freqs.size()]};
    }
    for (double rate : sweep_rates){
        std::cout << "Sweep rate = " << rate << " V/s" << std::endl;
        startTransient(V, n, p, Phi);
        double dt = fabs(VG_u - VG_l) / (n_VG - 1) / substeps / rate;
        for (int k = 0; k < n_VG; k++){
            int i = fromAccumulation(k);
            if (k > 0) moveGate(V, n, p, Phi, biasPoint(i) - phi_ms, substeps, dt);
            std::complex<double> C = smallSignalAC(x, V, n, p, column(TRANSIENT), 2 * pi * freqs.back(), V_ac, n_ac, p_ac);
            rate_data[i].push_back(C.real());
        }
        labels.push_back(label("C (sweep at ", rate, " V/s)"));
    }
    for (int i = 0; i < n_VG; i++){
        rate_data[i].push_back(vg_c_data[i].back());
    }
    labels.push_back("C deep depletion");
    save_table("csv_files/VG_C_sweep.csv", "VG", labels, rate_data);
}

// Transient 2: C-t. The gate is stepped (in a very short time) into inversion and then held, while
// the inversion layer builds up by generation and the capacitance recovers.
static void stepResponse(Matrix& V, Matrix& n, Matrix& p, Matrix& Phi){
    std::vector<std::complex<double>> V_ac, n_ac, p_ac;
    int tr = column(TRANSIENT);
    double VG_inv = p_type ? VG_u : VG_l;   // inversion end of the gate-voltage range

    // Times after the step grow geometrically
    std::vector<std::vector<double>> t_data;
    for (double t = t_first; t < t_last * growth; t *= growth) t_data.push_back({t});

    std::vector<std::string> labels;
    for (double fraction : step_fractions){
        double bias = fraction * VG_inv;
        std::cout << "Gate stepped to " << bias << " V" << std::endl;
        startTransient(V, n, p, Phi);
        int steps = std::max(1, (int)ceil(fabs(bias - phi_ms - V(0, tr)) / max_step));
        moveGate(V, n, p, Phi, bias - phi_ms, steps, t_ramp);
        double t = 0;
        for (auto& row : t_data){
            moveGate(V, n, p, Phi, bias - phi_ms, 1, row[0] - t);
            t = row[0];
            std::complex<double> C = smallSignalAC(x, V, n, p, tr, 2 * pi * freqs.back(), V_ac, n_ac, p_ac);
            row.push_back(C.real());
        }
        labels.push_back(label("C (VG stepped to ", bias, " V)"));
    }
    save_table("csv_files/t_C.csv", "t", labels, t_data);
}

int main(){
    // One column of V, n, p for each condition (see main.hpp)
    Matrix V(x_points + 1, n_columns), n(x_points + 1, n_columns), p(x_points + 1, n_columns);
    Matrix Phi(x_points + 1, 2);    // quasi-Fermi potentials of electrons (column 0) and holes (column 1), for the transient
    setupDevice(n, p);

    std::cout << "EQUILIBRIUM\n\n";
    rampGate(x, V, n, p, EQUILIBRIUM, 0 - phi_ms); // solve for equilibrium (VG = 0)
    save_x_based("csv_files/eq.csv", x, V, n, p, column(EQUILIBRIUM));

    std::cout << "BIAS SWEEP\n\n";
    std::vector<std::vector<double>> vg_c_data, vg_g_data;
    copyState(V, n, p, column(EQUILIBRIUM), column(DC));   // start the DC solution from equilibrium
    biasSweep(V, n, p, vg_c_data, vg_g_data);

    std::cout << "DEEP DEPLETION\n\n";
    deepDepletionSweep(V, n, p, vg_c_data);

    std::vector<std::string> c_labels, g_labels;
    for (double f : freqs){
        c_labels.push_back(label("C (f = ", f, " Hz)"));
        g_labels.push_back(label("G/w (f = ", f, " Hz)"));
    }
    c_labels.push_back("C deep depletion");
    save_table("csv_files/VG_C.csv", "VG", c_labels, vg_c_data);
    save_table("csv_files/VG_G.csv", "VG", g_labels, vg_g_data);

    std::cout << "FREQUENCY SWEEP\n\n";
    frequencySweep(V, n, p);

    std::cout << "TRANSIENT\n\n";
    sweepRateCV(V, n, p, Phi, vg_c_data);
    stepResponse(V, n, p, Phi);
    return 0;
}
