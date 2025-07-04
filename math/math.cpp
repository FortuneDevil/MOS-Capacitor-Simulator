#include "math.hpp"

// Bernoulli function = x/(e^x-1)
double B(double x){
    if(fabs(x) < 1e-6) return 1/(1+x/2);
    return x/(exp(x)-1);
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

        if (A(i, i) == 0){
                std::cout << "Can't be solved.\n";
                return x;
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

double max_relative_error(const Matrix& A, const Matrix& B) {
    if (A.size() != B.size()) {
        throw std::invalid_argument("Matrix sizes must match for error calculation");
    }

    double max_rel_error = 0.0;

    // Loop through each element of the matrices
    for (int i = 0; i < A.size(); ++i) {
        for (int j = 0; j < A.size(); ++j) {
            // Avoid division by zero in case of very small elements in B
            if (A(i, j) != 0) {
                double rel_error = std::fabs(1 - std::fabs(B(i, j)/A(i, j)));
                max_rel_error = std::max(max_rel_error, rel_error);  // Update the maximum relative error
            } else if (B(i, j) != 0){
                double rel_error = std::fabs(1 - std::fabs(A(i, j)/B(i, j)));
                max_rel_error = std::max(max_rel_error, rel_error);  // Update the maximum relative error
            }
        }
    }

    return max_rel_error;
}