#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include"h_files/matrix.h"
#include "h_files/nmdm.h"
#define n_VG 5//No. of points for VG
// Function to save data to CSV file
void save_to_csv(const char* filename, double* x, Matrix* V, Matrix* n, Matrix* p, int start_SC, int time) {
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Could not open file %s for writing.\n", filename);
        return;
    }

    // Write header
    fprintf(file, "x,V,n,p\n");

    // Write data
    for (int i = 0; i < x_points; i++) { // Loop through x array
        double voltage = (i < V->rows) ? V->a[i][time] : 0.0;
        double electron = (i >= start_SC && i < n->rows + start_SC) ? n->a[i - start_SC][time] : 0.0;
        double hole = (i >= start_SC && i < p->rows + start_SC) ? p->a[i - start_SC][time] : 0.0;
        fprintf(file, "%.6e,%.6e,%.6e,%.6e\n", x[i], voltage, electron, hole);
    }

    fclose(file);
    printf("Data saved to %s\n", filename);
}

int main(){
    double x[x_points+1]; //double epsilon[x_points+1];
    double t_ox =10e-9;
    double start_SC = 10;
    Matrix V = zero(x_points+1,2);
    //int start_SC = log(t_ox*(b_x-1)/c_x+1)/log(b_x)-1;
    //printf("%d\n",start_SC);
    Matrix n = zero(V.rows-start_SC-1, V.cols);
    Matrix p = zero(V.rows-start_SC-1, V.cols);
    Matrix V_eq = zero(V.rows,1);
    Matrix n_eq = zero(V_eq.rows-start_SC-1, 1);
    Matrix p_eq = zero(V_eq.rows-start_SC-1, 1);
    for (int i=0; i<=start_SC; i++){
		x[i] = i*t_ox/start_SC;
	}
	// spacing in the semiconductor is not uniform and is c_x b_x^i
	for (int i=start_SC+1; i<=x_points; i++){
		x[i] = t_ox + c_x *(pow(b_x,(i-start_SC))-1)/(b_x-1);
	}
    /*
    for(int i=0; i<=start_SC; i++){
        epsilon[i] = e_ox;
    }
    for(int i=start_SC+1; i<x_points+1; i++){
        epsilon[i] = e_si;
    }
    */
    n_eq.a[n_eq.rows-1][0] = pow(ni,2)/NA;
    p_eq.a[p_eq.rows-1][0] = NA;
    n_eq.a[0][0] = pow(ni,2)/NA;
    p_eq.a[0][0] = NA;

    consistent_eq(x,&V_eq,&n_eq,&p_eq,start_SC);
    //apply a dc bias
    double VG, VG_u=2, VG_l=-2; //VG_u -> upper bound of VG and VG_l for lower
    FILE *file = fopen("C_vs_VG.txt", "w");
    if (file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    fprintf(file, "VG\tC\n"); // Write header to file
    for(int i=0; i<n_VG; i++){
    VG = VG_l + (VG_u-VG_l)*i/n_VG;
    V=zero(V.rows,V.cols);
    V.a[i][0] = VG;
    for(int i=0; i<n.rows; i++){
        n.a[i][0] = n_eq.a[i][0];
        p.a[i][0] = p_eq.a[i][0];
    }
    consistent_dc(x,&V,&n,&p,&n_eq,&p_eq,start_SC);
    for(int i=0; i<V.rows; i++){
        V.a[i][1] = V.a[i][0];
    }
    for(int i=0; i<n.rows; i++){
        n.a[i][1] = n.a[i][0];
        p.a[i][1] = p.a[i][0];
    }
    V.a[0][1] += V.a[0][1]/100*(1e-3);
    consistent_ac(x, &V, &n, &p, &n_eq, &p_eq,start_SC);
    double C = e_ox/(x[1]-x[0])*(1-(V.a[1][1]-V.a[1][0])/(V.a[0][1]-V.a[0][0]))/1e2;
    fprintf(file, "%f\t%f\n", VG, C); // Write VG and C to file
    }
    fclose(file);
    printf("Data saved to C_vs_VG.txt\n");
    //VG = VAC+VDC
    //VAC = (VDC/100)sin(wt), if wt is small, then sin(wt)=wt
    // Save results to CSV file
    //save_to_csv("equilibrium.csv", x, &V_eq, &n_eq, &p_eq, start_SC,0);
    //save_to_csv("ac.csv", x, &V, &n, &p, start_SC,1);
    freeMatrix(&V);
    freeMatrix(&n);
    freeMatrix(&p);
    freeMatrix(&n_eq);
    freeMatrix(&p_eq);
    return 0;
}
