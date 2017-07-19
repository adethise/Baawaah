#include "cfd.h"



struct _problem* init_problem(){
  return malloc(sizeof(struct _problem));
}

void init_problem_physical(struct _problem* Problem, double CFL,double r_f, double L, double H, double Ls, double Hs, double h, double Q0, double tol, double nu,double Rey,double Str,double tau_max){
  (*Problem).CFL   = CFL;
  (*Problem).r_f     = r_f;
  (*Problem).L     = L;
  (*Problem).H     = H;
  (*Problem).Ls    = Ls;
  (*Problem).Hs    = Hs;
  (*Problem).h     = h;
  (*Problem).t     = 0;
  (*Problem).tau   = 0;
  (*Problem).tau_max = tau_max;
  (*Problem).Q0    = Q0;
  (*Problem).tol   = tol;
  (*Problem).nu    = nu;
  (*Problem).Rey   = Rey;
  (*Problem).Str   = Str;

  (*Problem).Um    = Q0/H;
  (*Problem).Uc    = 3.0/(2.0*(H-Hs)) * Q0;
  (*Problem).f     = Q0*Str/(H*H);
  (*Problem).Hc    = H-Hs;
  (*Problem).dt    = fmin(fabs(CFL*h/(*Problem).Uc),r_f*h*h/nu);
  (*Problem).dtau  = fabs((*Problem).dt*(*Problem).Um/L);
  (*Problem).tmax  = fabs((*Problem).tau_max*L/(*Problem).Um);
}


void init_problem_numerical(struct _problem* Problem, double phi,double t_snapshot,int flag_os,int flag_pres){
  (*Problem).Nx         = (*Problem).L    /(*Problem).h;
  (*Problem).Ny         = (*Problem).H    /(*Problem).h;
  (*Problem).Nx_p       = (*Problem).Nx +2;
  (*Problem).Ny_p       = (*Problem).Ny +2;
  (*Problem).Ntime      = (*Problem).tmax /(*Problem).dt;
  (*Problem).phi        = phi;
  (*Problem).t_snapshot = t_snapshot;
  (*Problem).flag_os    = flag_os;
  (*Problem).flag_pres  = flag_pres;
}

void init_problem_map(struct _problem* Problem){
  (*Problem).imap   = calloc((*Problem).Nx,sizeof(int));
  (*Problem).NHs = (*Problem).Hs/(*Problem).h;
  (*Problem).NLs = (*Problem).Ls/(*Problem).h;
  for(int i = 0 ; i < (*Problem).Nx ; i++){
      if( i < (*Problem).NLs ) (*Problem).imap[i] = (*Problem).NHs-1;
      else (*Problem).imap[i] = 0;
  }
}

void init_problem_vector_domain(struct _problem* Problem){
  (*Problem).omega       = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).psi         = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).u           = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).v           = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).f_old       = (double**) malloc((*Problem).Nx*sizeof(double*));
  (*Problem).R_res       = (double**) malloc((*Problem).Nx*sizeof(double*));

  (*Problem).P           = (double**) malloc(((*Problem).Nx+1)*sizeof(double*));
  (*Problem).R_res_pres  = (double**) malloc(((*Problem).Nx+1)*sizeof(double*));
  (*Problem).u_stag      = (double**) malloc((*Problem).Nx_p*sizeof(double*));
  (*Problem).v_stag      = (double**) malloc((*Problem).Nx_p*sizeof(double*));

  (*Problem).psi_in      = (double*) calloc((*Problem).Ny,sizeof(double));
  (*Problem).psi_out     = (double*) calloc((*Problem).Ny,sizeof(double));

  (*Problem).Re_h        = (double*) calloc((*Problem).Ntime+1,sizeof(double));
  (*Problem).Re_h_omega  = (double*) calloc((*Problem).Ntime+1,sizeof(double));
  (*Problem).Beta_CFL    = (double*) calloc((*Problem).Ntime+1,sizeof(double));


  for(int i=0 ; i < (*Problem).Nx ; i++){
    (*Problem).omega[i]      = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).psi[i]        = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).u[i]          = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).v[i]          = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).f_old[i]      = (double*) calloc((*Problem).Ny,sizeof(double));
    (*Problem).R_res[i]      = (double*) calloc((*Problem).Ny,sizeof(double));
  }
  for(int i=0 ; i < (*Problem).Nx+1 ; i++){
    (*Problem).P[i]          = (double*) calloc((*Problem).Ny+1,sizeof(double));
    (*Problem).R_res_pres[i] = (double*) calloc((*Problem).Ny+1,sizeof(double));
  }
  for(int i=0 ; i < (*Problem).Nx_p ; i++){
    (*Problem).u_stag[i]     = (double*) calloc((*Problem).Ny_p,sizeof(double));
    (*Problem).v_stag[i]     = (double*) calloc((*Problem).Ny_p,sizeof(double));
  }
}

void init_problem_poiseuille(struct _problem* Problem){
  for(int i = 0; i < (*Problem).Nx;i++ ){
    for(int j = (*Problem).NHs; j < (*Problem).Ny; j++ ){
      (*Problem).u[i][j]     = scalar_u_v_poiseuille(Problem,(j - (*Problem).NHs+1)*(*Problem).h);
    }
  }
    for(int i = 0; i < (*Problem).Nx;i++ ){
    for(int j = (*Problem).NHs; j < (*Problem).Ny; j++ ){
      (*Problem).omega[i][j]     = scalar_u_v_poiseuille_dy (Problem,(j - (*Problem).NHs+1)*(*Problem).h);
      //(*Problem).psi[i][j]       = scalar_u_v_poiseuille_int(Problem,(j - (*Problem).NHs+1)*(*Problem).h);
      // Maybe starting from a existing psi field can accelerate the convergence, but no time to investigate
    }
  }
}


void free_problem_vector_domain(struct _problem* Problem){
  for(int i=0 ; i < (*Problem).Nx ; i++){
    free((*Problem).omega[i]);
    free((*Problem).psi[i]);
    free((*Problem).u[i]);
    free((*Problem).v[i]);
    free((*Problem).f_old[i]);
    free((*Problem).R_res[i]);
  }
  for(int i=0 ; i < (*Problem).Nx_p ; i++){
    free((*Problem).u_stag[i]);
    free((*Problem).v_stag[i]);
  }
  for(int i=0 ; i < (*Problem).Nx_p-1 ; i++){
    free((*Problem).P[i]);
    free((*Problem).R_res_pres[i]);
  }

  free((*Problem).omega);
  free((*Problem).psi);
  free((*Problem).u);
  free((*Problem).v);

  free((*Problem).f_old);
  free((*Problem).R_res);

  free((*Problem).P);
  free((*Problem).u_stag);
  free((*Problem).v_stag);
  free((*Problem).R_res_pres);
  free((*Problem).psi_in);
  free((*Problem).psi_out);

  free((*Problem).imap);
  free((*Problem).Re_h);
  free((*Problem).Re_h_omega);
  free((*Problem).Beta_CFL);

  free(Problem);
}
