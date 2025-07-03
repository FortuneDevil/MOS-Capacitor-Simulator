#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"

Matrix createMatrix(int rows, int cols) {
    double **a = (double **)malloc(rows * sizeof(double *));
    for (int k = 0; k < rows; k++) {
        a[k] = (double *)malloc(cols * sizeof(double));
    }
    Matrix A;
    A.a = a;
    A.rows = rows;
    A.cols = cols;
    return A;
}

Matrix zero(int rows, int cols){
	Matrix A = createMatrix(rows,cols);
	for(int i=0; i<rows; i++){
		for(int j=0; j<cols; j++){
			A.a[i][j]=0;
		}
	}
	return A;
}

Matrix copy(Matrix* A){
	Matrix B = createMatrix(A->rows, A->cols);
	for(int i=0; i < A->rows; i++){
		for(int j=0; j < A->cols; j++){
			B.a[i][j] = A->a[i][j];
		}
	}
	return B;
}

Matrix extractrow(int row_num, Matrix* A) {
    Matrix B = createMatrix(1, A->cols);
    for (int i = 0; i < A->cols; i++) {
        B.a[0][i] = A->a[row_num][i];
    }
    return B;
}

Matrix extractcol(int col_num, Matrix* A) {
    Matrix B = createMatrix(A->rows, 1);
    for (int i = 0; i < A->rows; i++) {
        B.a[i][0] = A->a[i][col_num];
    }
    return B;
}

double extractelement(int i, int j, Matrix* A) {
    return A->a[i][j];
}

double add(double a, double b) {
    return a + b;
}

double sub(double a, double b) {
    return a - b;
}

Matrix addsubMatrix(Matrix* A, Matrix* B, double (*operation)(double a, double b)) {
    if (A->rows == B->rows && A->cols == B->cols) {
        Matrix C = createMatrix(A->rows, A->cols);
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->cols; j++) {
                C.a[i][j] = operation(A->a[i][j], B->a[i][j]);
            }
        }
        return C;
    } else {
        printf("Cannot add/subtract the given matrices.");
        return createMatrix(0, 0);
    }
}

Matrix scaleMatrix(double scalar, Matrix* A) {
    Matrix B = createMatrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            B.a[i][j] = scalar * A->a[i][j];
        }
    }
    return B;
}

void freeMatrix(Matrix* A) {
    for (int i = 0; i < A->rows; i++) {
        free(A->a[i]);
    }
    free(A->a);
}

void printMatrix(Matrix* A) {
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            printf("%15.4e ", A->a[i][j]);
        }
        printf("\n");
    }
}

Matrix transpose(Matrix* A) {
    Matrix B = createMatrix(A->cols, A->rows);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            B.a[j][i] = A->a[i][j];
        }
    }
    return B;
}

double det(Matrix* A) {
    if (A->rows == A->cols) {
        if (A->rows == 1) return A->a[0][0];
        if (A->rows == 2) return A->a[0][0] * A->a[1][1] - A->a[0][1] * A->a[1][0];
        double determinant = 0;
        for (int j = 0; j < A->cols; j++) {
            Matrix minor = minor_matrix(0, j, A);
            determinant += (j % 2 == 0 ? 1 : -1) * A->a[0][j] * det(&minor);
            freeMatrix(&minor);
        }
        return determinant;
    } else {
        printf("Not possible to calculate the determinant\n");
        return -1;
    }
}

Matrix minor_matrix(int i, int j, Matrix* A) {
    Matrix B = createMatrix(A->rows - 1, A->cols - 1);
    for (int row = 0; row < A->rows; row++) {
        for (int col = 0; col < A->cols; col++) {
            if (row != i && col != j) {
                int newRow = row < i ? row : row - 1;
                int newCol = col < j ? col : col - 1;
                B.a[newRow][newCol] = A->a[row][col];
            }
        }
    }
    return B;
}

Matrix cofactor(Matrix* A) {
    Matrix Ct = createMatrix(A->rows, A->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            Matrix minor = minor_matrix(i, j, A);
            Ct.a[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * det(&minor);
            freeMatrix(&minor);
        }
    }
    Matrix C = transpose(&Ct);
    freeMatrix(&Ct);
    return C;
}

Matrix inverse_nxn(Matrix* A){
	double Mat[A->rows][A->cols];
	double Factor,Factor1;
	Matrix InvMat=Identity(A->rows);
	for (int i = 0; i < A->rows; i++) {
		for (int j = 0; j < A->cols; j++) {
			Mat[i][j] = A->a[i][j];
		}
	}
	
	// Loop over all the rows
	for (int i=0;i<A->rows;i++){
		//Make diagonal elements non-zero by adding some other row.
		if(Mat[i][i]==0){
			int j=i;
			while(Mat[i][i]==0 && j<A->rows){
				Mat[i][i] += Mat[j][i];
				if(Mat[i][i]!=0){
					for(int k=0;k<A->cols;k++){
						if(k!=i){
						Mat[i][k]+=Mat[j][k];
						}
						InvMat.a[i][k]+=InvMat.a[j][k];
					}
				}
				j++;
			}
		}
		// Printing the matrix after making it an identity matrix
		//if still the factor1 is 0, this means no element in that column is non-zero, so det=0.
		if(Mat[i][i]==0){
		printf("It is singular matrix, inverse can't be calculated\n");
		return createMatrix(0,0);
		}
		Factor1 = Mat[i][i];
		for (int k=0;k<A->rows;k++){
			Mat[i][k] = Mat[i][k]/Factor1; // Make the diagonal elements to be one
			InvMat.a[i][k] = InvMat.a[i][k]/Factor1; // Use the same factor for inverse matrix
		}
		for (int j=0;j<A->rows;j++){ // Loop over all other rows
			if (i == j) {
				continue;
			}
			Factor = Mat[j][i];
			for (int k=0;k<A->cols;k++){ // Loop over columns
				Mat[j][k] = Mat[j][k]-Factor*Mat[i][k]; // Make the off-diagonal elements to be zero.
				InvMat.a[j][k] = InvMat.a[j][k]-Factor*InvMat.a[i][k];
			}
		}
	}
	return InvMat;
}


Matrix MultiplyMatrix(Matrix* A, Matrix* B) {
    if (A->cols != B->rows) {
        printf("Cannot multiply matrices: columns of A don't match rows of B\n");
        return createMatrix(0, 0);
    }
    Matrix C = createMatrix(A->rows, B->cols);
    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->cols; j++) {
            C.a[i][j] = 0;
            for (int k = 0; k < A->cols; k++) {
                C.a[i][j] += A->a[i][k] * B->a[k][j];
            }
        }
    }
    return C;
}

double trace(Matrix* A) {
    double sum = 0;
    for (int i = 0; i < A->rows && i < A->cols; i++) {
        sum += A->a[i][i];
    }
    if (A->rows != A->cols) printf("Warning: Trace is not defined, the trace value calculated might be incorrect\n");
    return sum;
}

double innerproduct(Matrix* A, Matrix* B) {
    Matrix Bt = transpose(B);
    Matrix product = MultiplyMatrix(&Bt, A);
    double result = trace(&product);
    freeMatrix(&Bt);
    freeMatrix(&product);
    return result;
}

void solution_A_nxnX_B(Matrix* A, Matrix* B, double* b) {
    if (B->cols == 1) {
    	Matrix inv_A = inverse_nxn(A);
        Matrix X = MultiplyMatrix(&inv_A, B);
        for (int i = 0; i < A->rows; i++) {
            b[i] = X.a[i][0];
        }
        freeMatrix(&X);
    } else {
        printf("Cannot calculate solution as B isn't a vector.\n");
    }
}

Matrix Identity(int rows) {
    Matrix I = createMatrix(rows, rows);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < rows; j++) {
            I.a[i][j] = (i == j) ? 1 : 0;
        }
    }
    return I;
}

double delta(int i, int j) {
    return (i == j) ? 1 : 0;
}

Matrix right_inverse(Matrix* A) {
    if (A->rows < A->cols) {
        Matrix At = transpose(A);
        Matrix inverse = MultiplyMatrix(A, &At);
        freeMatrix(&At);
        return inverse;
    } else {
        printf("Right inverse doesn't exist.\n");
        return createMatrix(0, 0);
    }
}

Matrix left_inverse(Matrix* A) {
    if (A->rows > A->cols) {
        Matrix At = transpose(A);
        Matrix inverse = MultiplyMatrix(&At, A);
        freeMatrix(&At);
        return inverse;
    } else {
        printf("Left inverse doesn't exist.\n");
        return createMatrix(0, 0);
    }
}

// Function to write the matrix to a file
void printMatrix_to_file(Matrix *m, FILE *file) {
    // Write the matrix elements to the file
    for (int i = 0; i < m->rows; i++) {
        for (int j = 0; j < m->cols; j++) {
            fprintf(file, "%15.4e ", m->a[i][j]);
        }
        fprintf(file, "\n");  // New line after each row
    }
}
