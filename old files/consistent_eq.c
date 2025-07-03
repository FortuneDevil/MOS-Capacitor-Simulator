#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
void consistent_eq(double* x,Matrix* V,Matrix* n,Matrix* p, int start_SC){
	int iter=0;
	Matrix temp_V, temp_n, temp_p;
	while(iter < max_iter){
		temp_V = copy(V);
		temp_n = copy(n);
		temp_p = copy(p);
		poisson(x, V, n, p, 0, start_SC);
		carrier_eq(x, V, n, p, start_SC);
		if(error(V,&temp_V,0) < tol_V  && error(n,&temp_n,0) < tol_n && error(p,&temp_p,0) < tol_p){
		//printf("%e,%e,%e\n",error(V,&temp_V,0),error(n,&temp_n,0), error(p,&temp_p,0));
			printf("Converged at iter=%d\n",iter);
			//printMatrix(n);
			//printMatrix(p);
			break;
		}
		iter++;
	}
	freeMatrix(&temp_V);
	freeMatrix(&temp_n);
	freeMatrix(&temp_p);
}