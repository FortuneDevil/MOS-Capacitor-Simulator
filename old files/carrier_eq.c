#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
double B(double x){
    if(fabs(x)<1e-6){
        return 1/(1+x/2);
    }
    return x/(exp(x)-1);
}

//Equilibrium carrier conc.
void carrier_eq(double* x,Matrix* V,Matrix* n,Matrix* p, int start_SC){
    //J i-1/2 = J i+1/2
    Matrix An = zero(n->rows,n->rows);
    Matrix Ap = zero(p->rows,p->rows);
    Matrix delta_n = zero(n->rows,1);
    Matrix delta_p = zero(p->rows,1);
    Matrix Fn = zero(n->rows,1);
    Matrix Fp = zero(p->rows,1);
    An.a[0][0] = 1; Ap.a[0][0] = 1;
    An.a[n->rows-1][n->rows-1] = 1; Ap.a[p->rows-1][p->rows-1] = 1;
    for (int i=start_SC+2; i<V->rows-1; i++) {
        int iSC = i-start_SC-1; //starts from 1, ends at V->rows-start_SC-3
        double delta_x_i = x[i+1]-x[i];
        double delta_x_i_minus = x[i]-x[i-1];
        An.a[iSC][iSC-1] = B((V->a[i-1][0]-V->a[i][0])/Vt)/delta_x_i_minus;
        An.a[iSC][iSC] = -(B((V->a[i][0]-V->a[i-1][0])/Vt)/delta_x_i_minus
                           + B((V->a[i][0]-V->a[i+1][0])/Vt)/delta_x_i );
        An.a[iSC][iSC+1] = B((V->a[i+1][0]-V->a[i][0])/Vt)/delta_x_i;
        Fn.a[iSC][0] = -(An.a[iSC][iSC-1]*n->a[iSC-1][0]
                        +An.a[iSC][iSC]*n->a[iSC][0]
                        +An.a[iSC][iSC+1]*n->a[iSC+1][0]);
        Ap.a[iSC][iSC-1] = B((V->a[i][0]-V->a[i-1][0])/Vt)/delta_x_i_minus;
        Ap.a[iSC][iSC] = -( B((V->a[i-1][0]-V->a[i][0])/Vt)/delta_x_i_minus 
                            + B((V->a[i+1][0]-V->a[i][0])/Vt)/delta_x_i );
        Ap.a[iSC][iSC+1] = B((V->a[i][0]-V->a[i+1][0])/Vt)/delta_x_i;
        Fp.a[iSC][0] = -(Ap.a[iSC][iSC-1]*p->a[iSC-1][0]
                        +Ap.a[iSC][iSC]*p->a[iSC][0]
                        +Ap.a[iSC][iSC+1]*p->a[iSC+1][0]);
    }
    /*
    printf("An:\n");
    printMatrix(&An);
    printf("n:\n");
    printMatrix(n);
    printf("p:\n");
    printMatrix(p);
    */
    gauss(&An,&delta_n,&Fn);
    *n = addsubMatrix(n,&delta_n,add);
    gauss(&Ap,&delta_p,&Fp);
    *p = addsubMatrix(p,&delta_p,add);
    freeMatrix(&An);
    freeMatrix(&delta_n);
    freeMatrix(&Fn);
    freeMatrix(&Ap);
    freeMatrix(&delta_p);
    freeMatrix(&Fp);
}
