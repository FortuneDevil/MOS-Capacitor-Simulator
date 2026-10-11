//parameters for delta_x[i]=c*b_x^i
#define b_x 1.01
#define c_x 1e-10

#define tol_V 1e-6 // Convergence tolerance
#define tol_n 1e-3 // Convergence tolerance
#define tol_p 1e-3 // Convergence tolerance
#define x_points 500 //no. of points in x
//#define t_points 0 //for no. of points in t=t_points.
#define max_iter 10000  // Maximum number of iterations to get self consistent soln
#define q 1.6e-19 
#define ni 1.5e10
#define kb 1.38e-23
#define T 300
#define ND 0
#define NA 1e17
#define Nc 3e19
#define Nv 1e19
#define un 1350
#define up 480
#define tau_n 1e-6
#define tau_p 1e-6
#define e_0 8.85e-14
#define e_si (11.7*e_0)
#define e_ox (3.9*e_0)
#define ni 1.5e10
#define Vt (kb*T/q)
#define freq 1e3 //omega
#define delta_t (1e-3/freq)
void carrier_ac(double* x,Matrix* V,Matrix* n,Matrix* p, Matrix* n_eq, Matrix* p_eq, int start_SC);
void carrier_dc(double* x,Matrix* V,Matrix* n,Matrix* p, Matrix* n_eq, Matrix* p_eq, int start_SC);
double B(double x);
void carrier_eq(double* x,Matrix* V,Matrix* n,Matrix* p, int start_SC);
void consistent_ac(double* x,Matrix* V,Matrix* n,Matrix* p,Matrix* n_eq,Matrix* p_eq, int start_SC);
void consistent_dc(double* x,Matrix* V,Matrix* n,Matrix* p,Matrix* n_eq,Matrix* p_eq, int start_SC);
void consistent_eq(double* x,Matrix* V,Matrix* n,Matrix* p, int start_SC);
double error(Matrix* A, Matrix* B, int col);
void gauss(Matrix* A, Matrix* x, Matrix* b);
void poisson(double* x,Matrix* V,Matrix* n,Matrix* p, int time, int start_SC);
