#include "poisson.hpp"

void poissonSolver(const std::vector<double>& x, Matrix& V, const Matrix& n, const Matrix& p, Condition cond){
    int N = V.size(); //  V is a column matrix with V.rows = N
    Matrix J(N, N);
    Matrix delta_v(N, 1);
    Matrix F(N, 1);
    int time;
    switch (cond) {
        case EQUILIBRIUM:
            time = 0;
            break;
        case DC:
            time = 1;
            break;
        case AC:
            time = 2;
            break;
        default:
            break;
    }
    J(0, 0) = 1.0;//starting point 
    double delta_x_ox = x[1] - x[0];  //spacing in oxide(const)
    double delta_x_ox_sq = pow(delta_x_ox, 2);//spacing in semiconductor

        for (int i = 1; i < start_SC; i++) {
        J(i, i - 1) = e_ox / (2*delta_x_ox_sq); //oxide i-1 index for the J matrix 
        J(i, i) = - ( e_ox / delta_x_ox_sq ); //oxide i index for the J matrix 
        J(i, i + 1) = e_ox / (2*delta_x_ox_sq); //oxide i+1 index for the J matrix 

        F(i, 0) = -(J(i, i - 1) * V(i - 1, time) // F matrix for oxide layer 
                  + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time));
    }
    // startSC point is assumed to be in boundary 
    double delta_x_sc = x[start_SC + 1] - x[start_SC]; //delta_x_sc is the difference between the first and second point of semiconductor

    J(start_SC, start_SC - 1) = e_ox / delta_x_ox; // initial conditions (borrowed from previous code)
    J(start_SC, start_SC) = -(e_ox / delta_x_ox + e_si / delta_x_sc);
    J(start_SC, start_SC + 1) = e_si / delta_x_sc;
    F(start_SC,0) = -(J(start_SC, start_SC - 1) * V(start_SC - 1, time)  // point in oxise so no  n , p terms
                  + J(start_SC, start_SC) * V(start_SC, time)
                  + J(start_SC, start_SC + 1) * V(start_SC + 1, time));

     for (int i = start_SC + 1; i < N - 1; i++) {
        double delta_x_i = x[i + 1] - x[i]; //distance between  x(i+1)-x(i)  
        double delta_x_i_sq = pow(delta_x_i, 2);
        double factor = e_si * sqrt(b_x) / delta_x_i_sq;//constant 
        //J matrix (trigonal one) and q is non zero
        J(i, i - 1) = factor * b_x;
        J(i, i) = - factor * (1 + b_x);
        J(i, i + 1) = factor;
        //F matrix  J.del = F
        // F = Vmatix + q;
         F(i, 0) = -(J(i, i - 1) * V(i - 1, time)
                    + J(i, i) * V(i, time)
                  + J(i, i + 1) * V(i + 1, time)
                  + q * (ND - NA + p(i, time) - n(i, time)));
     }
   // last condition where potential is zero at bulk 
      J(N - 1, N - 1) = 1.0;
      F(N - 1, 0) = 0.0;
   // gauss(J, delta_v, F);
   delta_v = solver(J, F);

    for (int i = 0; i < N; i++) {
        V(i, time) += delta_v(i, 0);
    }
}