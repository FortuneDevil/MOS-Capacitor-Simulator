#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
double error(Matrix* A, Matrix* B, int col){
    double x=0,max=0;
    for(int i=0; i<A->rows; i++){
        if(A->a[i][col] != 0){
            x = (A->a[i][col]-B->a[i][col])/A->a[i][col];
        } else if(B->a[i][col] != 0){
            x = (A->a[i][col]-B->a[i][col])/B->a[i][col];
        } else x = 0;
        
        if(x > max){
            max = x;
        }
    }
    return max;
}