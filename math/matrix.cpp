#include "matrix.hpp"

// Constructor: initialize with rows and cols, filled with zeros
Matrix::Matrix(int r, int c): data(r * c, 0.0), rows(r), cols(c){}

// Copy constructor
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols), data(other.data) {}

// Identity Matrix
Matrix Matrix::identity(int size){
    Matrix I(size, size);
    for(int i = 0; i < size; i++){
        I(i, i) = 1.0;
    }
    return I;
}

int Matrix::size(){
    return rows;
}

int Matrix::size() const {
    return rows;
}


void Matrix::swapRows(int row1, int row2) {
    for (int col = 0; col < cols; ++col) {
        std::swap(data[row1 * cols + col], data[row2 * cols + col]);
    }
}


// Assignment operator
Matrix& Matrix::operator=(const Matrix& B) {
    if(this == &B)
        return *this;
    
    rows = B.rows;
    cols = B.cols;
    data = B.data;

    return *this;
}

// Access operator (modifiable)
double& Matrix::operator()(int i, int j){
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw std::out_of_range("Index out of range");
    return data[i * cols + j];
}

// Access operator (const version)
const double& Matrix::operator()(int i, int j) const{
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw std::out_of_range("Index out of range");
    return data[i * cols + j];
}

// Matrix Addition
Matrix Matrix::operator+(const Matrix& B) const{
    if (rows != B.rows || cols != B.cols)
        throw std::invalid_argument("Matrix sizes didn't match for addition");
    Matrix C(rows, cols);
    for(int i = 0; i < rows * cols; i++){
        C.data[i] = data[i] + B.data[i];
    }
    return C;
}

// Matrix Subtraction
Matrix Matrix::operator-(const Matrix& B) const{
    if (rows != B.rows || cols != B.cols)
        throw std::invalid_argument("Matrix sizes didn't match for subtraction");
    Matrix C(rows, cols);
    for(int i = 0; i < rows * cols; i++){
        C.data[i] = data[i] - B.data[i];
    }
    return C;
}

// Matrix * scalar
Matrix Matrix::operator*(double scalar) const{
    Matrix C(rows, cols);
    for (int i = 0; i < rows * cols; i++)
        C.data[i] = data[i] * scalar;
    return C;
}

// Matrix * Matrix
Matrix Matrix::operator*(const Matrix& B) const{
    if (cols != B.rows)
        throw std::invalid_argument("Matrix dimensions didn't match for multiplication");
    
    Matrix C(rows, B.cols);

    for(int i = 0; i < rows; i++){
        for (int j = 0; j < B.cols; j++){
            double sum = 0;
            for (int k = 0; k < cols; k++){
                sum += (*this)(i, k) * B(k, j);
            }
            C(i, j) = sum;
        }
    }
    return C;
}

// Transpose of a Matrix
Matrix Matrix::transpose() const{
    Matrix AT(cols, rows);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            AT(j, i) = (*this)(i, j);
        }
    }
    return AT;
}

// Printing the Matrix
void Matrix::printMatrix(std::ostream& out) const{
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            out << std::setw(12) << std::setprecision(5) << std::scientific << (*this)(i, j) << " ";
        }
        out << "\n";
    }
}


std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
    matrix.printMatrix(os);  // Call the printMatrix method to print the matrix
    return os;
}

// scalar * matrix
Matrix operator*(double scalar, const Matrix& M){
    return M * scalar;
}