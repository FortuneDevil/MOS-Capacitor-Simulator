#include "matrix.hpp"

int main() {
    Matrix A = Matrix::identity(3);
    std::cout << "Printing A:\n";
    A.printMatrix();
    Matrix B(3, 3);

    B(0, 0) = 2.0;
    B(1, 1) = 3.0;
    B(2, 2) = 4.0;
    B(2, 1) = 1e-9;
    
    std::cout << "Printing B:\n";
    B.printMatrix();

    Matrix C = A + B;
    std::cout << "Printing C:\n";
    C.printMatrix();

    Matrix D = C * 2.0;
    std::cout << "Printing D:\n";
    D.printMatrix();

    Matrix E = 2.0 * C;
    std::cout << "Printing E:\n";
    E.printMatrix();

    Matrix F = B * C;
    std::cout << "Printing F:\n";
    F.printMatrix();

    Matrix G = F; // copies F to G
    std::cout << "G(=F) : \n" << G;

    std::cout << "Printing BT:\n";
    B.transpose().printMatrix();

    std::cout << "This also prints matrix A:\n" << A;
    return 0;
}
