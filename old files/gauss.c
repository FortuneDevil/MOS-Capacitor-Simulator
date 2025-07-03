#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
void gauss(Matrix* A, Matrix* x, Matrix* b){
	double n = A->rows;
    //Forward Elimination
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
            double factor = A->a[j][i] / A->a[i][i];
            for (int k = i; k < n; k++) {
                A->a[j][k] -= A->a[i][k] * factor;
            }
            b->a[j][0] -= b->a[i][0] * factor;
        }
	}
    // Backward Substitution
    for (int i = n - 1; i >= 0; i--) {
        x->a[i][0] = b->a[i][0];
        for (int j = i + 1; j < n; j++) {
            x->a[i][0] -= A->a[i][j] * x->a[j][0];
        }
        x->a[i][0] /= A->a[i][i];
    }
}