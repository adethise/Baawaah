#include "cfd.h"

double test_Qfunc_const(struct _problem* Problem,double t){
  return 1.0;
}
void test_omega_domainFill(struct _problem* Problem){
  for(int i = 0; i < (*Problem).Nx ; i++) {
    for(int j = 0; j < (*Problem).imax_map[i] ; j++ ) (*Problem).omega[i][j] = 1;
  }
}
void test_psi_boundaryFill(struct _problem* Problem, double (*Q)(double) ){
// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).psi[i][0] = Q((*Problem).t);
// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ){
    if   (i < (*Problem).NLs ){ (*Problem).psi[i][(*Problem).NHs-1] = 2; }
    else                        (*Problem).psi[i][(*Problem).Ny -1] = 2;
  }
// Down boundary Side
  for(int i = (*Problem).NHs; i < (*Problem).Ny; i++ ) (*Problem).psi[(*Problem).NLs-1][i] = 2;
// Inflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).NHs-1 ; i++) (*Problem).psi[0][i]               = (*Problem).psi[1][i];
// Outflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).Ny -1 ; i++) (*Problem).psi[(*Problem).Nx-1][i] = (*Problem).psi[(*Problem).Nx-2][i];
}

void test_omega_boundaryFill(struct _problem* Problem){
// Upper boundary
  for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][0] = 2;
// Down boundary - Left - Right
  for(int i = 0; i < (*Problem).Nx ; i++ ){
    if   (i < (*Problem).NLs ){ (*Problem).omega[i][(*Problem).NHs-1] = 2; }
    else                        (*Problem).omega[i][(*Problem).Ny -1] = 2;
  }
// Down boundary Side
  for(int i = (*Problem).NHs; i < (*Problem).Ny; i++ ) (*Problem).omega[(*Problem).NLs-1][i] = 2;
// Corner boundary
  (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] = 3;
  (*Problem).omega[(*Problem).NLs-1][(*Problem).Ny -1] = 3;
// Inflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).NHs-1 ; i++) (*Problem).omega[0][i]             = 4;
// Outflow Boundary - Natural Condition
  for(int i = 1; i < (*Problem).Ny -1 ; i++) (*Problem).omega[(*Problem).Nx-1][i] = 4;
}
