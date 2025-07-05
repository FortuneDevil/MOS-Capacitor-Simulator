#include <cmath>
#include "matrix.hpp"
using namespace std;
//constants i will put later

void poisson(const std::vector<double>& x, Matrix& V, const Matrix& n, const Matrix& p, int time, int start_SC) {
    int N = V.size(); // assuming V is a column matrix with V.rows = N
    Matrix J(N, N);
    Matrix delta_v(N, 1);
    Matrix F(N, 1);

    J(0, 0) = 1.0;
    double delta_x_ox = x[1] - x[0];  //spacing in oxide(const)
    double delta_x_ox_sq = pow(delta_x_ox, 2);//spacing in semiconductor
   // double delta_prev = ;
    //double delta_next = ;

    for (int i = 1; i < start_SC; i++) {
        J(i, i - 1) = e_ox / (2*delta_x_ox_sq);
        J(i, i) = - ( e_ox / delta_x_ox_sq );
        J(i, i + 1) = e_ox / (2*delta_x_ox_sq);
      // delta_prev = ;
      // delta_next = ;
      //J(i,i-1) = 1/((delta_next-delta_prev)(delta_next));
      //J(i,i) = -({1/(delta_next-delta_prev)}{1/(delta_prev)+1/(delta_next)});
      //J(i,i+1) = 

        F(i, 0) = -(J(i, i - 1) * V(i - 1, time)
                  + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time));
    }

    double delta_x_sc = x[start_SC + 1] - x[start_SC];
    double deno = (delta_x_ox + delta_x_sc) / 2.0;
    J(start_SC, start_SC - 1) = e_ox / delta_x_ox / deno;
    J(start_SC, start_SC) = -(e_ox / delta_x_ox + e_si / delta_x_sc) / deno;
    J(start_SC, start_SC + 1) = e_si / delta_x_sc / deno;

    for (int i = start_SC + 1; i < N - 1; i++) {
        int iSC = i - start_SC - 1;
        double delta_x_i = x[i + 1] - x[i];
        double delta_x_i_minus = x[i] - x[i - 1];
        double deno_sc = (delta_x_i + delta_x_i_minus);

        J(i, i - 1) =   e_si  /( delta_x_i_minus*deno_sc);
        J(i, i) = - ( e_si / (deno_sc))(1/delta_x_i+1/delta_x_i_minus);
                + q * (n(iSC, time) + p(iSC, time)) / Vt;
        J(i, i + 1) =  e_si / (delta_x_i*deno_sc);

        F(i, 0) = -(J(i, i - 1) * V(i - 1, time)
                    + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time)
                  + q * (ND - NA + p(iSC, time) - n(iSC, time)));
    }

    //    for (int i = start_SC + 1; i < N - 1; i++) {
    //    int iSC = i - start_SC - 1;
    //    double delta_x_i = x[i + 1] - x[i];
    //    double delta_x_i_minus = x[i] - x[i - 1];
    //    double deno_sc = delta_x_i * delta_x_i_minus * (delta_x_i + delta_x_i_minus);

    //    J(i, i - 1) = 2 * e_si * delta_x_i / deno_sc;
     //   J(i, i) = -2 * e_si * (delta_x_i + delta_x_i_minus) / deno_sc
     //           + q * (n(iSC, time) + p(iSC, time)) / Vt;
     //   J(i, i + 1) = 2 * e_si * delta_x_i_minus / deno_sc;

    //    F(i, 0) = -(J(i, i - 1) * V(i - 1, time)
    //              - 2 * e_si * (delta_x_i + delta_x_i_minus) / deno_sc * V(i, time)
    //              + J(i, i + 1) * V(i + 1, time)
    //              + q * (ND - NA + p(iSC, time) - n(iSC, time)));
    //}

    J(N - 1, N - 1) = 1.0;

    gauss(J, delta_v, F);

    for (int i = 0; i < N; i++) {
        V(i, time) += delta_v(i, 0);
    }
}
