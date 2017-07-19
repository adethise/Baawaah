/*
 *  CFD - Linux Version of the Main
 *  by Son Tran
 */

#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <unistd.h>
#include "cfd.h"
#include <getopt.h>


void deadstop_exit(struct _problem* Problem){
  fprintf(stderr, "[DEADSTOP EXIT] \n");
  fprintf(stderr, "case diagnostic DEADSTOP  \n");
  fprintf(stderr, "    Binary Value for each problem\n");
  fprintf(stderr, "       Mesh Reynold u,v   - 1 \n");
  fprintf(stderr, "       Mesh Reynold omega - 2 \n");
  fprintf(stderr, "       Beta - CFL         - 4 \n");
  fprintf(stderr, "Emergency Exit Procedure Initiate\n");
  fprintf(stderr, "Problem occured at tau = %d\n",(int) ( (*Problem).tau * 100 ));
  fprintf(stderr, "Saving current data under  _%d_\n",(int) ( (*Problem).tau * 100 ));
  print_problem_data    (Problem);
  print_problem_pressure(Problem);
  //print_problem_diag    (Problem); // Need to be investigated
  free_problem_vector_domain(Problem);
  exit(-1);
}

void print_problem_diag(struct _problem* Problem){
    FILE* file_diag = fopen("data/CFD_DIAG.txt","w");
    if(file_diag == NULL){ fprintf(stderr,"File error\n"); exit(1);}
    for(int i = 0; i < (*Problem).Ntime+1 ; i++ ){
    fprintf(file_diag,"%f %5.16f %5.16f %5.16f \n",i*(*Problem).dtau,(*Problem).Re_h[i],(*Problem).Re_h_omega[i],(*Problem).Beta_CFL[i]);
    }
    fclose(file_diag);
}

void print_problem_pressure(struct _problem* Problem){
  char buff_P_name[50];
  char buff_R_name[50];
  sprintf(buff_P_name ,"data/CFD_P_%d.txt",(int) ( (*Problem).tau * 100 ));
  sprintf(buff_R_name ,"data/CFD_R_pres_%d.txt",(int) ( (*Problem).tau * 100 ));
  FILE* file_P = fopen(buff_P_name,"w");
  FILE* file_R = fopen(buff_R_name,"w");
  if(file_P == NULL || file_R == NULL){ fprintf(stderr,"File error\n"); exit(1);}
    for(int j = 0 ; j < (*Problem).Ny_p-1 ; j++){
      for(int i = 0; i < (*Problem).Nx_p-1 ; i++){
        fprintf(file_P,"%5.16f "   ,(*Problem).P[i][j]);
        fprintf(file_R,"%5.16f "   ,(*Problem).R_res_pres[i][j]);
    }
    fprintf(file_P ,"\n");
    fprintf(file_R ,"\n");
  }
  fclose(file_P);
  fclose(file_R);
}

void print_problem_data(struct _problem* Problem){
  char buff_omega_name[50];
  char buff_psi_name[50];
  char buff_u_name[50];
  char buff_v_name[50];
  char buff_R_name[50];
  sprintf(buff_omega_name ,"data/CFD_omega_%d.txt",(int) ( (*Problem).tau * 100 ));
  sprintf(buff_psi_name   ,"data/CFD_psi_%d.txt"  ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_u_name     ,"data/CFD_u_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_v_name     ,"data/CFD_v_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  sprintf(buff_R_name     ,"data/CFD_R_%d.txt"    ,(int) ( (*Problem).tau * 100 ));
  FILE* file_omega = fopen(buff_omega_name,"w");
  FILE* file_psi   = fopen(buff_psi_name  ,"w");
  FILE* file_u     = fopen(buff_u_name    ,"w");
  FILE* file_v     = fopen(buff_v_name    ,"w");
  FILE* file_R     = fopen(buff_R_name    ,"w");

  if(file_omega == NULL || file_psi == NULL || file_u == NULL || file_v == NULL){ fprintf(stderr,"File error\n"); exit(1);}
    for(int j = 0 ; j < (*Problem).Ny ; j++){
      for(int i = 0; i < (*Problem).Nx ; i++){
        fprintf(file_omega,"%5.16f "   ,(*Problem).omega[i][j]);
        fprintf(file_psi  ,"%5.16f "   ,(*Problem).psi  [i][j]);
        fprintf(file_u    ,"%5.16f "   ,(*Problem).u    [i][j]);
        fprintf(file_v    ,"%5.16f "   ,(*Problem).v    [i][j]);
        fprintf(file_R    ,"%5.16f "   ,(*Problem).R_res[i][j]);
    }
    fprintf(file_omega ,"\n");
    fprintf(file_psi   ,"\n");
    fprintf(file_u     ,"\n");
    fprintf(file_v     ,"\n");
    fprintf(file_R     ,"\n");
  }
  fclose(file_omega);
  fclose(file_psi);
  fclose(file_u);
  fclose(file_v);
  fclose(file_R);
}

int main(int argc,char* argv[]){
  // Flow specification
  double nu    =   1e-6;

  // Numerical parameter
  double CFL        =   0.2 ;
  double r_f        =   0.5 ;
  double tau_max    =   1.0 ;
  double Rey        =   100 ;
  double Str        =   0.5 ;
  double h          =   0.02;
  double t_snapshot =   0.25;
  int    flag_os    =   0   ;
  int    flag_pres  =   0   ;

  // Domain parameter
  double H     =     1.0 ;
  double L     =   15.0*H;
  double Ls    =   5.0*H ;
  double Hs    =   H/2.0 ;
  double Q0    =   Rey*nu;
  double sign  =   1.0;
  // Computation parameter
  //
  char option;
  while ((option = getopt(argc, argv, "boph:R:C:r:t:s:")) != EOF) {
          switch (option) {
          case 'b':
              sign = -1.0;
              break;
          case 'o':
              flag_os = 1;
              break;
          case 'p':
              flag_pres = 1;
              break;
          case 'h':
              sscanf(optarg, "%lf", &h);
              break;
          case 'R':
              sscanf(optarg, "%lf", &Rey);
              break;
          case 'C':
              sscanf(optarg, "%lf", &CFL);
              break;
          case 'r':
              sscanf(optarg, "%lf", &r_f);
              break;
          case 't':
              sscanf(optarg, "%lf", &tau_max);
              break;
          case 's':
              sscanf(optarg, "%lf", &t_snapshot);
              break;
          }
  }
  if( h <= 0 || h >= 1 || t_snapshot > 1 || tau_max < 0 || Rey <= 0 || CFL <= 0){printf("[DEADSTOP] Input Error\n");exit(-1);}

  //
  Q0    =  sign*Rey*nu;
  double phi   =   1.98;
  double tol   =   1e-4;

  // ########
  fprintf(stderr, "CFD - Oscillating flow in a channel with abrupt step\n");
  fprintf(stderr, "by S. Tran - J. Demey \n");

  struct _problem* Problem = init_problem();
  init_problem_physical(Problem,CFL,r_f,L,H,Ls,Hs,h,Q0,tol,nu,Rey,Str,tau_max);
  init_problem_numerical(Problem,phi,t_snapshot,flag_os,flag_pres);
  init_problem_map(Problem);
  init_problem_vector_domain(Problem);
  init_problem_poiseuille(Problem);
  fprintf(stderr, "CFL: %f Fourier: %f Tau: %f Reynold: %f Strouhal: %f\n",(*Problem).CFL,(*Problem).r_f,(*Problem).tau_max,(*Problem).Rey,(*Problem).Str);
  fprintf(stderr, "Average Velocity: %f Q0: %f Grid size h: %f Timestep: %f Frequency: %f Snapshot: %f\n",(*Problem).Um,(*Problem).Q0,(*Problem).h,(*Problem).dt,(*Problem).f,(*Problem).t_snapshot);
  fprintf(stderr, "NX: %d NY: %d NLs: %d NHs: %d NTime: %d \n",(*Problem).Nx,(*Problem).Ny,(*Problem).NLs,(*Problem).NHs,(*Problem).Ntime);
  fprintf(stderr,"|Option| Facing Step Mode: %d  | Oscillating Mode: %d  | Pressure Mode: %d |\n",Q0<0,(*Problem).flag_os,(*Problem).flag_pres);
  printf("Simulation Starting\n");


  boundary_psi_update(Problem,functionQ);
  poisson_inner_psi_iterator(Problem);
  boundary_omega_update(Problem);
  inner_u_v_compute(Problem);


  print_problem_data(Problem);
  print_problem_pressure(Problem);


  // ---Code Benchmarking-------
  struct timespec start, finish;
  double elapsed;
  clock_gettime(CLOCK_MONOTONIC, &start);
  // ---------------------------

  integration_omega(Problem);

  // ---Code Benchmarking-------
  clock_gettime(CLOCK_MONOTONIC, &finish);
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  // ---------------------------
  if((*Problem).flag_pres == 1){
  u_v_stag(Problem);
  boundary_pression_ghost_update(Problem);
  boundary_pression_ghost_in_out(Problem);
  poisson_inner_pres_iterator(Problem);
  }

  print_problem_data(Problem);
  print_problem_diag(Problem);
  print_problem_pressure(Problem);
  printf("Simulation Done\n");


  printf("Time Elapsed: %f s\n",elapsed);
  printf("Final File save under _%d_\n",(int) ( (*Problem).tau * 100 ));
  free_problem_vector_domain(Problem);
  return 0;
}
