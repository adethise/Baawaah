/*
 * CFD.h
 * CFD General Structure and methods titles
 */
 #include <stdlib.h>
 #include <math.h>

   //Domain Physical Propiety
   struct _problem {
   double CFL;
   double r_f;
   double L;
   double H;
   double Ls;
   double Hs;
   double Hc;
   double h;
   double dt;
   double dtau;
   double t;
   double tau;
   double tmax;
   double tau_max;
   double Q0;
   double nu;
   double Um;
   double Uc;
   double Rey;
   double Str;
   double f;
   //Domain Numerical property
   int    Nx;
   int    Ny;
   int    Nx_p;
   int    Ny_p;
   int    NLs;
   int    NHs;
   int    Ntime;
   double phi;
   double e_max;
   double tol;
   int   * imap;
   double* Re_h;
   double* Re_h_omega;
   double* Beta_CFL;

   int    flag_os;
   int    flag_pres;

   double t_snapshot;
   //Domain Data of size N_x * N_y
   double **omega;
   double **psi;
   double **u;
   double **v;
   double **f_old;
   double **R_res;
   // Staggered Part - Or The quite hard and long part
   double **u_stag;
   double **v_stag;
   double *psi_in;
   double *psi_out;
   double **P;
   double **R_res_pres;

 };
void             print_problem_diag           (struct _problem* Problem);
void             print_problem_pressure       (struct _problem* Problem);
void             print_problem_data           (struct _problem* Problem);
void             deadstop_exit                (struct _problem* Problem);
/* ###################################
 *  CFD Main
 * ###################################
 */
struct _problem* init_problem();
void             init_problem_physical        (struct _problem* Problem, double CFL, double r_f, double L, double H, double Ls, double Hs, double h, double Q0, double tol, double nu,double Rey,double Str, double tau_max);
void             init_problem_numerical       (struct _problem* Problem, double phi, double t_snapshot,int flag_os,int flag_pres);
void             init_problem_map             (struct _problem* Problem);
void             init_problem_vector_domain   (struct _problem* Problem);
void             init_problem_poiseuille      (struct _problem* Problem);
void             free_problem_vector_domain   (struct _problem* Problem);
/* ###################################
 *  CFD Numerical
 * ###################################
 */
void             scalar_rhs                    (struct _problem* Problem, int i, int j, double du);
double           scalar_psi_star_compute       (struct _problem* Problem,int i,int j);
double           scalar_psi_compute            (struct _problem* Problem,int i,int j);
double           scalar_psi_r_compute          (struct _problem* Problem,int i, int j);
double           scalar_pres_star_compute      (struct _problem* Problem,int i,int j);
double           scalar_Ax                     (struct _problem* Problem, int i, int j);
double           scalar_Ay                     (struct _problem* Problem, int i, int j);
double           scalar_u_v_poiseuille         (struct _problem* Problem,double y);
double           scalar_u_v_poiseuille_dy      (struct _problem* Problem,double y);
double           scalar_u_v_poiseuille_int     (struct _problem* Problem,double y);
void             inner_u_v_compute             (struct _problem* Problem);
void             u_v_stag                      (struct _problem* Problem);
void             boundary_psi_update           (struct _problem* Problem, double (*Q)(struct _problem*) );
void             boundary_omega_update         (struct _problem* Problem);
void             boundary_omega_dwdx_update    (struct _problem* Problem);
void             boundary_pression_ghost_update(struct _problem* Problem);
void             boundary_pression_ghost_in_out(struct _problem* Problem);
void             boundary_pression_ghost_corner(struct _problem* Problem,int flag);
void             inner_psi_update              (struct _problem* Problem);
double           inner_psi_error_compute       (struct _problem* Problem);
void             inner_pres_update             (struct _problem* Problem);
double           inner_pres_error_compute      (struct _problem* Problem);
void             poisson_inner_psi_iterator    (struct _problem* Problem);
void             poisson_inner_pres_iterator   (struct _problem* Problem);
int              diagnose_check                (struct _problem* Problem,int i,int j, int ktime);
void             snapshot                      (struct _problem* Problem,int k);
/* ###################################
 *  CFD Integrator
 * ###################################
 */
double           functionQ                    (struct _problem* Problem);
void             first_iteration_omega        (struct _problem* Problem);
void             first_time_integration       (struct _problem* Problem);
void             integration_omega            (struct _problem* Problem);
