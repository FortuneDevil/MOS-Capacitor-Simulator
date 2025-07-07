#include "carrier.hpp"

void carrierDC(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p){
    int N = V.size() - start_SC - 1;
    Matrix F_n(N, 1), F_p(N, 1);
    Matrix delta_n(N, 1), delta_p(N, 1);
    Matrix Jac_n(N, N), Jac_p(N, N);
    std::vector<double> phi (N, 0);
    int correctionFactor = start_SC + 1;
    for (int i = start_SC + 1; i < V.size(); i++){
        phi[i - correctionFactor] = V(i, 1) / Vt;
    }
    double factor1n = sqrt(b_x) * Vt * mu_n;
    double factor1p = sqrt(b_x) * Vt * mu_p;
    for (int i = 1; i < N - 1; i++){
        double delta_x_i = (x[i + 1] - x[i]);
        double delta_x_i_minus = (x[i] - x[i - 1]);
        double nfactor = factor1n / delta_x_i;
        double pfactor = factor1p / delta_x_i;
        int gi = i + correctionFactor;  // global index


        Jac_n(i, i + 1) = - nfactor * B(phi[i + 1] - phi[i]) / delta_x_i;
        Jac_n(i, i) = delta_x_i / tau_n + nfactor * (B(phi[i] - phi[i - 1]) / delta_x_i_minus 
                            + B(phi[i] - phi[i + 1]) / delta_x_i);
        Jac_n(i, i - 1) = - nfactor * B(phi[i - 1] - phi[i]) / delta_x_i_minus;

        double constant = - delta_x_i * n(i, 0) / tau_n;

        F_n(i, 0) = -(Jac_n(i, i + 1) * n(gi + 1, 1) + Jac_n(i, i) * n(gi, 1) 
                    + Jac_n(i, i - 1) * n(gi - 1, 1) + constant);


        // Now do the same for p
        Jac_p(i, i + 1) = - pfactor * B(phi[i] - phi[i + 1]) / delta_x_i;
        Jac_p(i, i) = delta_x_i / tau_p + pfactor * (B(phi[i - 1] - phi[i]) / delta_x_i_minus 
                            + B(phi[i + 1] - phi[i]) / delta_x_i);
        Jac_p(i, i - 1) = - pfactor * B(phi[i] - phi[i + 1]) / delta_x_i_minus;

        constant = - delta_x_i * p(i, 0) / tau_p;

        F_p(i, 0) = -(Jac_p(i, i + 1) * p(gi + 1, 1) + Jac_p(i, i) * p(gi, 1) 
                    + Jac_p(i, i - 1) * p(gi - 1, 1) + constant);
    }

    // Boundary conditions
    F_n(N - 1, 0) = 0; F_p(N - 1, 0) = 0;
    Jac_n(N - 1, N - 1) = 1; Jac_p(N - 1, N - 1) = 1;
    Jac_n(0, 0) = - exp(phi[1] - phi[0]);
    Jac_n(0, 1) = 1; 
    Jac_p(0, 1) = -1; 
    Jac_p(0, 0) = exp(phi[1] - phi[0]);

    // Solve for n, p
    delta_n = solver(Jac_n, F_n);
    delta_p = solver(Jac_p, F_p);

    for(int i = 0; i < N; i++){
        int gi = i + correctionFactor;  // global index
        n(gi, 1) += delta_n(i, 0);
        p(gi, 1) += delta_p(i, 0);
    }
}