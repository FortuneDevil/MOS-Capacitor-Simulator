#include "math.hpp"

// Bernoulli function = x/(e^x-1)
double B(double x){
    return bernoulli(x);
}

// Derivative of the Bernoulli function
double dB(double x){
    return dBernoulli(x);
}

// Solve for x in Ax = b using Gaussian Elimination
Matrix solver(Matrix &A, Matrix &b){
    int n = A.size();
    Matrix x(n, 1);

    std::vector<double> temp_row(n);    // Temporary row for swapping
    
    // Forward Elimination
    for (int i = 0; i < n; i++){

        // Find row with largest element in column i for pivoting
        int max_row = i;
        for (int j = i + 1; j < n; j++){
            if (fabs(A(j, i)) > fabs(A(max_row, i))){
                max_row = j;
            }
        }

        // Swap if max_row is not this row
        if (max_row != i) {
            // Swap rows in A
            A.swapRows(i, max_row);

            // Swap corresponding entries in b
            b.swapRows(i, max_row);
        }

        /*
        if (A(i, i) == 0){
                std::cout << "Can't be solved due to row "<< i << std::endl;
                A(i, i) = 1e-10;
        } 
        */
        if (std::abs(A(i, i)) < 1e-100) {
		std::cerr << "WARNING: Pivot nearly zero at row " << i << ", A(i,i)=" << A(i,i) << std::endl;
		A(i, i) = 1e-100;
	}
        for (int j = i + 1; j < n; j++){
            double factor = A(j, i)/A(i, i);
            for (int k = i; k < n; k++){
                A(j, k) -= A(i, k) * factor;
            }
            b(j, 0) -= b(i, 0) * factor;
        }
    }
    // Backward Substitution
    for (int i = n - 1; i >= 0; i--){
        x(i, 0) = b(i, 0);
        for (int j = i + 1; j < n; j++){
            x(i, 0) -= A(i, j) * x(j, 0);
        }
        x(i, 0) /= A(i,i);
    }
    return x;
}

std::vector<double> solveFluxSystem(const std::vector<double>& up, const std::vector<double>& down,
                                    const std::vector<double>& sink, const std::vector<double>& source, double last_value){
    int n = sink.size();
    std::vector<double> y(n), e(n), beta(n);

    // Forward Elimination. After it, node k obeys: up[k] * y[k+1] - (down[k] + e[k]) * y[k] = -beta[k]
    e[0] = sink[0];
    beta[0] = source[0];
    for (int k = 1; k < n - 1; k++){
        double pivot = down[k - 1] + e[k - 1];
        e[k] = sink[k] + up[k - 1] * e[k - 1] / pivot;
        beta[k] = source[k] + down[k - 1] * beta[k - 1] / pivot;
    }

    // Backward Substitution
    y[n - 1] = last_value;
    for (int k = n - 2; k >= 0; k--){
        y[k] = (up[k] * y[k + 1] + beta[k]) / (down[k] + e[k]);
    }
    return y;
}

double max_relative_error(const Matrix& A, const Matrix& B, int col) {
    if (A.size() != B.size()) {
        throw std::invalid_argument("Matrix sizes must match for error calculation");
    }

    double max_rel_error = 0.0;

    // Loop through each element of the matrices
    for (int i = 0; i < A.size(); ++i) {
        // Avoid division by zero in case of very small elements in B
        if (A(i, col) != 0) {
            double rel_error = std::fabs(1 - std::fabs(B(i, col)/A(i, col)));
            max_rel_error = std::max(max_rel_error, rel_error);  // Update the maximum relative error
        } else if (B(i, col) != 0){
            double rel_error = std::fabs(1 - std::fabs(A(i, col)/B(i, col)));
            max_rel_error = std::max(max_rel_error, rel_error);  // Update the maximum relative error
        }
    }
    return max_rel_error;
}
