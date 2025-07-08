#include "carrier.hpp"

void carrierEq(const std::vector<double>& x, const Matrix& V, Matrix& n, Matrix& p){
    int N = V.size() - start_SC - 1;
    Matrix F_n(N, 1), F_p(N, 1);
    Matrix delta_n(N, 1), delta_p(N, 1);
    Matrix Jac_n(N, N), Jac_p(N, N);
    std::vector<double> phi(N,0);
    for (int i = start_SC + 1; i < V.size(); i++){
        phi[i - start_SC -1] = V(i, 0) / Vt; 
    }

   //Boundary conditions
    F_n(N - 1, 0) = 0; F_p(N - 1, 0) = 0;
    Jac_n(N - 1, N - 1) = 1; Jac_p(N - 1, N - 1) = 1;
    Jac_n(0, 0) = - exp(phi[1] - phi[0]);
    Jac_n(0, 1) = 1; 
    Jac_p(0, 1) = -1; 
    Jac_p(0, 0) = exp(phi[1] - phi[0]);

   for(int  i = 1 ; i < N-1 ; i++ ){
   Jac_n(i, i + 1) = B(phi[ i + 1] - phi[ i]);
   Jac_n(i, i) = -(B(phi[ i]-phi[ i + 1]) + b_x * B(phi[i]-phi[i - 1]));
   Jac_n(i,i - 1) = b_x*B(phi[i - 1]-phi[i]); 

   F_n(i, 0) = Jac_n(i, i + 1) * n(start_SC + 1 + i + 1, 0) + Jac_n(i, i) * n(start_SC + 1 + i,0) + Jac_n(i, i - 1) * n(start_SC + 1 + i - 1, 0); 

   Jac_p(i, i + 1) = b_x*B(phi[i] - phi[i + 1]);
   Jac_p(i, i) = -(B(phi[i - 1]-phi[i]) + b_x*B(phi[i + 1]-phi[i]));
   Jac_p(i,i - 1) = B(phi[i]-phi[i - 1]); 

   F_p(i, 0) = Jac_p(i, i + 1) * p(start_SC + 1 + i + 1, 0) + Jac_p(i, i) * p(start_SC + 1 + i, 0) + Jac_p(i,  i - 1) * p(start_SC + 1 + i - 1, 0); 

   }
   
   //std::cout << "Carrier EQ: " << std::endl;
   delta_n = solver(Jac_n,F_n);
   delta_p = solver(Jac_p,F_p);

   for (int i = 0; i < N - 1; i++) {
        n(start_SC + 1 + i, 0) += dampingFactor * delta_n(i, 0);
        p(start_SC + 1 + i, 0) += dampingFactor * delta_p(i, 0);
    }
    
}
