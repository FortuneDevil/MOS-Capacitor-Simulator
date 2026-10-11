#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
void poisson(double* x,Matrix* V,Matrix* n,Matrix* p, int time, int start_SC){
    Matrix J = zero(V->rows,V->rows);
    Matrix delta_v = zero(V->rows,1);
    Matrix F = zero(V->rows,1);
    //J delta_x = F // this F is -F in our normal notations
    J.a[0][0] = 1;
    double delta_x_ox = x[1]-x[0];
	double delta_x_ox_sq = pow(delta_x_ox,2);
    for (int i=1; i<start_SC; i++) {
        J.a[i][i-1] = e_ox/delta_x_ox_sq;
		J.a[i][i] = -2*e_ox/delta_x_ox_sq;
		J.a[i][i+1] = e_ox/delta_x_ox_sq;
        F.a[i][0] = -(J.a[i][i-1]*V->a[i-1][time]
                    +J.a[i][i]*V->a[i][time]
                    +J.a[i][i+1]*V->a[i+1][time]);
    }
    //At interface, assuming start_SC in oxide.
    double delta_x_sc = x[start_SC+1]-x[start_SC];
    double deno = (delta_x_ox +delta_x_sc)/2;
    J.a[start_SC][start_SC-1] = e_ox/delta_x_ox/deno;
    J.a[start_SC][start_SC] = -(e_ox/delta_x_ox+e_si/delta_x_sc)/deno;
    J.a[start_SC][start_SC+1] = e_si/delta_x_sc/deno;
    //In semiconductor
    for (int i=start_SC+1; i<V->rows-1;i++){
        int iSC = i-start_SC-1;
        double delta_x_i = x[i+1]-x[i];
        double delta_x_i_minus = x[i]-x[i-1];
        double deno_sc=(delta_x_i*delta_x_i_minus*(delta_x_i+delta_x_i_minus));
        J.a[i][i-1] = 2*e_si*delta_x_i/deno_sc;
		J.a[i][i] = -2*e_si*(delta_x_i+delta_x_i_minus)/deno_sc
                    +q*(n->a[iSC][time]+p->a[iSC][time])/Vt;
		J.a[i][i+1] = 2*e_si*delta_x_i_minus/deno_sc;
        F.a[i][0] = -(J.a[i][i-1]*V->a[i-1][time]
                    - 2*e_si*(delta_x_i+delta_x_i_minus)/deno_sc*V->a[i][time]
                    + J.a[i][i+1]*V->a[i+1][time]
                    + q*(ND - NA + p->a[iSC][time] - n->a[iSC][time]));
    }
    J.a[V->rows-1][V->rows-1] = 1;
    gauss(&J,&delta_v,&F);
    for(int i=0; i<V->rows; i++){
        V->a[i][time] += delta_v.a[i][0];
    }
    freeMatrix(&J);
    freeMatrix(&delta_v);
    freeMatrix(&F);
}
