#include "math.hpp"

int main(){
// Example: 3x3 matrix
    Matrix A(3, 3);
    A(0, 0) = 3; A(0, 1) = 2; A(0, 2) = -1;
    A(1, 0) = 2; A(1, 1) = -1; A(1, 2) = 3;
    A(2, 0) = 1; A(2, 1) = 3; A(2, 2) = 2;

    // Right-hand side vector b
    Matrix b(3, 1);
    b(0, 0) = 1;
    b(1, 0) = 2;
    b(2, 0) = 3;

    // Solve Ax = b using Gaussian Elimination
    try {
        Matrix solution = solver(A, b);
        std::cout << "Solution: \n" << solution;
    } catch (const std::runtime_error& e) {
        std::cout << e.what() << std::endl;
    }

    return 0;
}