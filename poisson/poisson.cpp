#include <cmath>
#include "../math/math.hpp"
#include "../math/matrix.hpp"
#include "poisson.hpp"
//using namespace std;
//constants i will put later
//const double epsilon_0 = 8.854e-12;       // F/m
//const double kappa_ox = 3.9;
//const double kappa_si = 11.7;

//const double e_ox = epsilon_0 * kappa_ox; // 3.453e-11 F/m
//const double e_si = epsilon_0 * kappa_si; // 1.035e-10 F/m
//const double q = 1.9e-19;
//const double Vt = 0.5;

void poisson(const std::vector<double>& x, Matrix& V, const Matrix& n, const Matrix& p, int time, int start_SC) {
    int N = V.size(); //  V is a column matrix with V.rows = N
    Matrix J(N, N);
    Matrix delta_v(N, 1);
    Matrix F(N, 1);

    J(0, 0) = 1.0;//starting point 
    double delta_x_ox = x[1] - x[0];  //spacing in oxide(const)
    double delta_x_ox_sq = pow(delta_x_ox, 2);//spacing in semiconductor

        for (int i = 1; i < start_SC; i++) {
        J(i, i - 1) = e_ox / (2*delta_x_ox_sq); //oxide i-1 index for the J matrix 
        J(i, i) = - ( e_ox / delta_x_ox_sq ); //oxide i index for the J matrix 
        J(i, i + 1) = e_ox / (2*delta_x_ox_sq); //oxide i+1 index for the J matrix 
      // delta_prev = ;
      // delta_next = ;
      //J(i,i-1) = 1/((delta_next-delta_prev)(delta_next));
      //J(i,i) = -({1/(delta_next-delta_prev)}{1/(delta_prev)+1/(delta_next)});
      //J(i,i+1) = 

        F(i, 0) = -(J(i, i - 1) * V(i - 1, time) // F matrix for oxide layer 
                  + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time));
    }
    // startSC point is assumed to be in boundary 
    double delta_x_sc = x[start_SC + 1] - x[start_SC]; //delta_x_sc is the difference between the first and second point of semiconductor

    double deno = (delta_x_ox + delta_x_sc) / 2.0; // deno is the average of both the lengths
    J(start_SC, start_SC - 1) = e_ox / delta_x_ox / deno; // initial conditions (borrowed from previous code)
    J(start_SC, start_SC) = -(e_ox / delta_x_ox + e_si / delta_x_sc) / deno;
    J(start_SC, start_SC + 1) = e_si / delta_x_sc / deno;
    F(start_SC,0) = -(J(start_SC, start_SC - 1) * V(start_SC - 1, time)  // point in oxise so no  n , p terms
                  + J(start_SC, start_SC) * V(start_SC, time)
                  + J(start_SC, start_SC + 1) * V(start_SC + 1, time));

     for (int i = start_SC + 1; i < N - 1; i++) {
        double delta_x_i = x[i + 1] - x[i]; //distance between  x(i+1)-x(i)  
        double delta_x_i_minus = x[i] - x[i - 1];//distance between  x(i)-x(i-1) 
        double deno_sc = (delta_x_i + delta_x_i_minus);//constant 
        //J matrix (trigonal one) and q is non zero
         J(i, i - 1) =   e_si  /( delta_x_i_minus*deno_sc);
         J(i, i) = - ( e_si / (deno_sc))*(1/delta_x_i+1/delta_x_i_minus);
                + q * (n(i, time) + p(i, time)) / Vt;
         J(i, i + 1) =  e_si / (delta_x_i*deno_sc);
        //F matrix  J.del = F
        // F = Vmatix + q;
         F(i, 0) = -(J(i, i - 1) * V(i - 1, time)
                    + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time)
                  + q * (ND - NA + p(i, time) - n(i, time)));
     }
   // last condition where potential is zero at bulk 
      J(N - 1, N - 1) = 1.0;
      F(N - 1,0) = 0.0;
   // gauss(J, delta_v, F);
   delta_v = solver(J,F);

    for (int i = 0; i < N; i++) {
        V(i, time) += delta_v(i, 0);
    }

     


}