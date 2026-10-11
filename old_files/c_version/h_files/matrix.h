#ifndef MATRIX_H
#define MATRIX_H

// Structure for Matrix
typedef struct {
    double **a;  // Pointer to a 2D array of doubles
    int rows;    // Number of rows in the matrix
    int cols;    // Number of columns in the matrix
} Matrix;

// Function declarations
Matrix createMatrix(int rows, int cols);
Matrix zero(int rows, int cols);
Matrix copy(Matrix* V);
Matrix extractrow(int row_num, Matrix* A);
Matrix extractcol(int col_num, Matrix* A);
double extractelement(int i, int j, Matrix* A);
double add(double a, double b);
double sub(double a, double b);
Matrix addsubMatrix(Matrix* A, Matrix* B, double (*operation)(double a, double b));
Matrix scaleMatrix(double scalar, Matrix* A);
void freeMatrix(Matrix* A);
void printMatrix(Matrix* A);
Matrix transpose(Matrix* A);
double det(Matrix* A);
Matrix minor_matrix(int i, int j, Matrix* A);
Matrix cofactor(Matrix* A);
Matrix inverse_nxn(Matrix* A);
Matrix MultiplyMatrix(Matrix* A, Matrix* B);
double trace(Matrix* A);
double innerproduct(Matrix* A, Matrix* B);
void solution_A_nxnX_B(Matrix* A, Matrix* B, double* b);
Matrix Identity(int rows);
double delta(int i, int j);
Matrix right_inverse(Matrix* A);
Matrix left_inverse(Matrix* A);
// End function declarations

#endif // MATRIX_H