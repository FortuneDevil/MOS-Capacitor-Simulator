#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

class Matrix {
private:
    std::vector<double> data;
    int rows, cols;

public:
    Matrix(int rows, int cols);
    Matrix(const Matrix& other);
    
    static Matrix identity(int size);

    int size();
    int size() const;

    void swapRows(int row1, int row2);
    void copyColumn(const int a, const int b); // copies column a to b

    Matrix& operator=(const Matrix& B);
    double& operator()(int i, int j);
    const double& operator()(int i, int j) const;
    Matrix operator+(const Matrix& B) const;
    Matrix operator-(const Matrix& B) const;
    Matrix operator*(double scalar) const; // Scalar Multiplication
    Matrix operator*(const Matrix& B) const;    // Matrix Multiplication

    Matrix transpose() const;
    
    void printMatrix(std::ostream& out = std::cout) const;
    // Overload operator<< for matrix printing
    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix);
};

Matrix operator*(double scalar, const Matrix& M);
#endif