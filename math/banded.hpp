#ifndef BANDED_H
#define BANDED_H

#include <vector>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <algorithm>

// Square banded matrix: A(i, j) can be non-zero only for -kl <= j - i <= ku.
// kl extra super-diagonals are stored to hold the fill-in created by pivoting.
template <typename T>
class BandMatrix {
private:
    int n, kl, ku, width;
    std::vector<T> data;

public:
    BandMatrix(int n, int kl, int ku)
        : n(n), kl(kl), ku(ku), width(2 * kl + ku + 1), data((size_t)n * width, T(0)) {}

    int size() const { return n; }
    int lower() const { return kl; }
    int upper() const { return ku; }

    T& operator()(int i, int j){
        if (i < 0 || i >= n || j < 0 || j >= n || j - i < -kl || j - i > ku + kl)
            throw std::out_of_range("BandMatrix index outside the band");
        return data[(size_t)i * width + (j - i + kl)];
    }
};

// Solve A x = b by Gaussian elimination with scaled partial pivoting, working only inside the band.
// T can be a real or a complex type. A and b are modified.
template <typename T>
std::vector<T> solveBanded(BandMatrix<T>& A, std::vector<T>& b){
    int n = A.size(), kl = A.lower(), ku = A.upper();

    // Largest entry of each row, so that pivots are compared on rows of equal scale
    std::vector<long double> scale(n, 0);
    for (int i = 0; i < n; i++){
        for (int j = std::max(0, i - kl); j <= std::min(n - 1, i + ku); j++){
            scale[i] = std::max(scale[i], (long double)std::abs(A(i, j)));
        }
        if (scale[i] == 0) throw std::runtime_error("solveBanded: empty row in matrix");
    }

    // Forward Elimination
    for (int i = 0; i < n; i++){
        int last_row = std::min(n - 1, i + kl);
        int last_col = std::min(n - 1, i + ku + kl);

        int max_row = i;
        for (int r = i + 1; r <= last_row; r++){
            if (std::abs(A(r, i)) / scale[r] > std::abs(A(max_row, i)) / scale[max_row]){
                max_row = r;
            }
        }
        if (max_row != i){
            for (int j = i; j <= last_col; j++) std::swap(A(i, j), A(max_row, j));
            std::swap(b[i], b[max_row]);
            std::swap(scale[i], scale[max_row]);
        }
        if (std::abs(A(i, i)) == 0) throw std::runtime_error("solveBanded: singular matrix");

        for (int r = i + 1; r <= last_row; r++){
            T factor = A(r, i) / A(i, i);
            if (factor == T(0)) continue;
            for (int j = i; j <= last_col; j++){
                A(r, j) -= factor * A(i, j);
            }
            b[r] -= factor * b[i];
        }
    }

    // Backward Substitution
    std::vector<T> x(n);
    for (int i = n - 1; i >= 0; i--){
        T sum = b[i];
        for (int j = i + 1; j <= std::min(n - 1, i + ku + kl); j++){
            sum -= A(i, j) * x[j];
        }
        x[i] = sum / A(i, i);
    }
    return x;
}

#endif
