#ifndef BLOCK3_H
#define BLOCK3_H

#include <array>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>

// Small 3-vectors and 3x3 matrices, real or complex. They hold the three coupled unknowns of one node
// (potential, electron and hole quasi-Fermi potentials) and the three equations
// (Poisson, electron continuity, hole continuity).
template <typename T> using Vec3 = std::array<T, 3>;
template <typename T> using Mat3 = std::array<std::array<T, 3>, 3>;

template <typename T>
Vec3<T> zeroVec3(){
    Vec3<T> v;
    v.fill(T(0));
    return v;
}

template <typename T>
Mat3<T> zeroMat3(){
    Mat3<T> M;
    for (auto& row : M) row.fill(T(0));
    return M;
}

template <typename T>
Mat3<T> add3(const Mat3<T>& A, const Mat3<T>& B){
    Mat3<T> C;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) C[i][j] = A[i][j] + B[i][j];
    return C;
}

template <typename T>
Vec3<T> add3(const Vec3<T>& a, const Vec3<T>& b){
    return Vec3<T>{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

template <typename T>
Mat3<T> mul3(const Mat3<T>& A, const Mat3<T>& B){
    Mat3<T> C = zeroMat3<T>();
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++) C[i][j] += A[i][k] * B[k][j];
    return C;
}

template <typename T>
Vec3<T> mul3(const Mat3<T>& A, const Vec3<T>& v){
    Vec3<T> w;
    for (int i = 0; i < 3; i++) w[i] = A[i][0] * v[0] + A[i][1] * v[1] + A[i][2] * v[2];
    return w;
}

// Inverse by Gauss-Jordan elimination with scaled partial pivoting
template <typename T>
Mat3<T> inverse3(Mat3<T> A){
    Mat3<T> inv = zeroMat3<T>();
    long double scale[3];
    for (int i = 0; i < 3; i++){
        inv[i][i] = T(1);
        scale[i] = std::max({(long double)std::abs(A[i][0]), (long double)std::abs(A[i][1]), (long double)std::abs(A[i][2])});
    }
    for (int i = 0; i < 3; i++){
        int max_row = i;
        for (int r = i + 1; r < 3; r++){
            if (std::abs(A[r][i]) / scale[r] > std::abs(A[max_row][i]) / scale[max_row]) max_row = r;
        }
        std::swap(A[i], A[max_row]);
        std::swap(inv[i], inv[max_row]);
        std::swap(scale[i], scale[max_row]);

        T pivot = A[i][i];
        for (int j = 0; j < 3; j++){
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }
        for (int r = 0; r < 3; r++){
            if (r == i) continue;
            T factor = A[r][i];
            for (int j = 0; j < 3; j++){
                A[r][j] -= factor * A[i][j];
                inv[r][j] -= factor * inv[i][j];
            }
        }
    }
    return inv;
}

// Solves the block-tridiagonal system of three coupled conservation laws on nodes 0..N-1:
//   F(i) - F(i-1) - local[i] * X[i] = -source[i],   with fluxes F(i) = up[i] * X[i+1] - down[i] * X[i]
// between nodes i and i+1, no flux on the left of node 0, and X[N-1] = 0.
// This is the block version of solveFluxSystem: the elimination never subtracts the (large, nearly
// equal) fluxes from each other, so small local terms are not lost to round-off.
template <typename T>
std::vector<Vec3<T>> solveBlockFluxSystem(const std::vector<Mat3<T>>& up, const std::vector<Mat3<T>>& down,
                                          const std::vector<Mat3<T>>& local, const std::vector<Vec3<T>>& source){
    int N = local.size();
    std::vector<Mat3<T>> pivot_inv(N, zeroMat3<T>());
    std::vector<Vec3<T>> beta(N, zeroVec3<T>());

    // Forward Elimination. After it, node i obeys: up[i] * X[i+1] - (down[i] + E_i) * X[i] = -beta[i]
    Mat3<T> E = local[0];
    beta[0] = source[0];
    for (int i = 0; i < N - 1; i++){
        pivot_inv[i] = inverse3(add3(down[i], E));
        if (i + 1 < N - 1){
            E = add3(local[i + 1], mul3(mul3(E, pivot_inv[i]), up[i]));
            beta[i + 1] = add3(source[i + 1], mul3(down[i], mul3(pivot_inv[i], beta[i])));
        }
    }

    // Backward Substitution
    std::vector<Vec3<T>> X(N, zeroVec3<T>());
    for (int i = N - 2; i >= 0; i--){
        X[i] = mul3(pivot_inv[i], add3(mul3(up[i], X[i + 1]), beta[i]));
    }
    return X;
}

#endif
