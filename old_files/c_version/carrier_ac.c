#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
//DC carrier conc.
void carrier_ac(double* x,Matrix* V,Matrix* n,Matrix* p, Matrix* n_eq, Matrix* p_eq, int start_SC){
    Matrix An = zero(n->rows,n->rows);
    Matrix Ap = zero(p->rows,p->rows);
    Matrix Bn = zero(n->rows,1);
    Matrix Bp = zero(p->rows,1);
    Matrix delta_n = zero(n->rows,1);
    Matrix delta_p = zero(p->rows,1);
    An.a[0][0] = 1; An.a[1][0] = -1/exp((V->a[1][1]-V->a[0][1])/Vt);
    Ap.a[0][0] = 1; Ap.a[1][0] = -exp((V->a[1][1]-V->a[0][1])/Vt);
    An.a[n->rows-1][n->rows-1] = 1; Ap.a[p->rows-1][p->rows-1] = 1;
    for (int i=start_SC+2; i<V->rows-1; i++) {
        int iSC = i-start_SC-1; //starts from 1, ends at V->rows-start_SC-3
        double delta_x_i = x[i+1]-x[i];
        double delta_x_i_minus = x[i]-x[i-1];
        double delta_x_i_min_half = 1/sqrt(b_x) * delta_x_i;
        An.a[iSC][iSC-1] = -Vt*un/delta_x_i_min_half/2*B((V->a[i-1][1]-V->a[i][1])/Vt)/delta_x_i_minus;
        An.a[iSC][iSC] = (Vt*un/delta_x_i_min_half/2
                        *(B((V->a[i][1]-V->a[i-1][1])/Vt)/delta_x_i_minus
                         +B((V->a[i][1]-V->a[i+1][1])/Vt)/delta_x_i)+1/tau_n/2+1/delta_t);
        An.a[iSC][iSC+1] = -Vt*un/delta_x_i_min_half/2*B((V->a[i+1][1]-V->a[i][1])/Vt)/delta_x_i;
        Bn.a[iSC][0] = -(An.a[iSC][iSC-1]*n->a[iSC][iSC-1]
                        +An.a[iSC][iSC]*n->a[iSC][iSC]
                        +An.a[iSC][iSC+1]*n->a[iSC][iSC+1]
                        -n_eq->a[iSC][0]/tau_n/2 - n->a[iSC][0]/delta_t);
        Ap.a[iSC][iSC-1] = Vt*up/delta_x_i_min_half/2*B((V->a[i][1]-V->a[i-1][1])/Vt)/delta_x_i_minus;
        Ap.a[iSC][iSC] = -(Vt*up/delta_x_i_min_half/2
                         *(B((V->a[i-1][1]-V->a[i][1])/Vt)/delta_x_i_minus
                         +B((V->a[i+1][1]-V->a[i][1])/Vt)/delta_x_i)+1/tau_p/2+1/delta_t);
        Ap.a[iSC][iSC+1] = Vt*up/delta_x_i_min_half/2*B((V->a[i][1]-V->a[i+1][1])/Vt)/delta_x_i;
        Bp.a[iSC][0] = -(Ap.a[iSC][iSC-1]*p->a[iSC][iSC-1]
                        +Ap.a[iSC][iSC]*p->a[iSC][iSC]
                        +Ap.a[iSC][iSC+1]*p->a[iSC][iSC+1]
                        -p_eq->a[iSC][0]/tau_p/2 - p->a[iSC][0]/delta_t);
    }
    Bn.a[0][0] = (An.a[0][0]*n->a[0][1] - An.a[1][1]*n->a[1][1]);
    Bp.a[0][0] = -(Ap.a[0][0]*p->a[0][1] - Ap.a[1][0]*p->a[1][1]);
    Bn.a[n->rows-1][0] = n->a[n->rows-1][1];
    Bp.a[p->rows-1][0] = p->a[p->rows-1][1];
    gauss(&An,&delta_n,&Bn);
    gauss(&Ap,&delta_p,&Bp);
    for(int i=0; i<n->rows; i++){
        n->a[i][1] += delta_n.a[i][0];
        p->a[i][1] += delta_p.a[i][0];
    }
    freeMatrix(&An);
    freeMatrix(&Bn);
    freeMatrix(&Ap);
    freeMatrix(&Bp);
}

