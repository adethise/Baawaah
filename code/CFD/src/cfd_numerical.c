#include "cfd.h"
#include "stdio.h"
#include <stddef.h>

double scalar_psi_star_compute(struct _problem* Problem,int i,int j){
    return 0.25*(
                 (*Problem).omega[i][j] * ((*Problem).h*(*Problem).h)
                 + ((*Problem).psi[i+1][j] + (*Problem).psi[i-1][j])
                 + ((*Problem).psi[i][j+1] + (*Problem).psi[i][j-1])    );
}

double scalar_psi_compute(struct _problem* Problem,int i,int j){
    return    (1.0-(*Problem).phi) * (*Problem).psi[i][j]
    +      (*Problem).phi  * scalar_psi_star_compute(Problem,i,j);
}

double scalar_psi_r_compute(struct _problem* Problem,int i, int j){
    return   (*Problem).omega[i][j] +
    (  (*Problem).psi[i+1][j]   + (*Problem).psi[i-1][j]
     - 4*(*Problem).psi[i][j]
     +(*Problem).psi[i]  [j+1] + (*Problem).psi[i]  [j-1]
     )/pow((*Problem).h,2);
}
double scalar_pres_star_compute(struct _problem* Problem,int i,int j){
    return 0.25*(
                 (scalar_Ax(Problem,i+1,j) - scalar_Ax(Problem,i,j) - scalar_Ay(Problem,i,j+1) + scalar_Ay(Problem,i,j) ) * ((*Problem).h)
                 + ((*Problem).P[i+1][j] + (*Problem).P[i-1][j])
                 + ((*Problem).P[i][j+1] + (*Problem).P[i][j-1])    );
}

double scalar_pres_compute(struct _problem* Problem,int i,int j){
    return    (1.0-(*Problem).phi) * (*Problem).P[i][j]
    +      (*Problem).phi  * scalar_pres_star_compute(Problem,i,j);
}

double scalar_pres_r_compute(struct _problem* Problem,int i, int j){
    return  (scalar_Ax(Problem,i+1,j) - scalar_Ax(Problem,i,j) - scalar_Ay(Problem,i,j+1) + scalar_Ay(Problem,i,j)             )/(*Problem).h
    + ((*Problem).P[i+1][j] + (*Problem).P[i-1][j] - 4*(*Problem).P[i][j] + (*Problem).P[i][j+1] + (*Problem).P[i][j-1])/pow((*Problem).h,2);
}

double scalar_u_v_poiseuille     (struct _problem* Problem,double y){
    return 6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc) *  ( y            - (y*y)/(*Problem).Hc              );
}
double scalar_u_v_poiseuille_dy  (struct _problem* Problem,double y){
    return -6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc)  * ( 1.0          - (2.0*y)/(*Problem).Hc            );
}
double scalar_u_v_poiseuille_int (struct _problem* Problem,double y){
    return 6.0*functionQ(Problem)/((*Problem).Hc*(*Problem).Hc)  * ((1.0/2.0)*y*y - (1.0/3.0)*(y*y*y)/(*Problem).Hc) + functionQ(Problem);
}

void inner_u_v_compute(struct _problem* Problem){
    for(int i= 1 ; i < (*Problem).Nx-1; i++){
        for(int j=(*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
            // by centred finite differences
            (*Problem).u    [i][j] =     ((*Problem).psi[i][j+1]-(*Problem).psi[i][j-1])/(2.0*(*Problem).h);
            (*Problem).v    [i][j] =    -((*Problem).psi[i+1][j]-(*Problem).psi[i-1][j])/(2.0*(*Problem).h);
        }
    }
}

double scalar_Ax(struct _problem* Problem, int i, int j){
    double part1, part2;
    part1 = (*Problem).u_stag[i][j]*((*Problem).u_stag[i+1][j]-(*Problem).u_stag[i-1][j])/(2.0*(*Problem).h);
    part2 = ((*Problem).v_stag[i+1][j]+(*Problem).v_stag[i+1][j+1]+(*Problem).v_stag[i][j+1]+(*Problem).v_stag[i][j])/4;
    return part1+part2*((*Problem).u_stag[i][j+1]-(*Problem).u_stag[i][j-1])/(2.0*(*Problem).h);
}

double scalar_Ay(struct _problem* Problem, int i, int j){
    double part1, part2;
    part1 = (*Problem).v_stag[i][j]*((*Problem).v_stag[i][j+1]-(*Problem).v_stag[i][j-1])/(2.0*(*Problem).h);
    part2 = ((*Problem).u_stag[i+1][j]+(*Problem).u_stag[i+1][j+1]+(*Problem).u_stag[i][j+1]+(*Problem).u_stag[i][j])/4;
    return part1+part2*((*Problem).v_stag[i+1][j]-(*Problem).v_stag[i-1][j])/(2.0*(*Problem).h);
}

void u_v_stag(struct _problem* Problem){
    for(int i= 1 ; i < (*Problem).Nx_p-2; i++){
        for(int j=(*Problem).imap[i]; j < (*Problem).Ny_p-2; j++){
            (*Problem).u_stag    [i][j] =     ((*Problem).psi[i-1][j]-(*Problem).psi[i-1][j-1])/(*Problem).h;
            (*Problem).v_stag    [i][j] =    -((*Problem).psi[i][j-1]-(*Problem).psi[i-1][j-1])/(*Problem).h;
        }
    }
    // Upper
    for(int i = 1; i < (*Problem).Nx_p-2 ; i++ )
        (*Problem).u_stag[i][(*Problem).Ny_p-1] = -0.2*(15*(*Problem).u_stag[i][(*Problem).Ny_p-2] -5.0*(*Problem).u_stag[i][(*Problem).Ny_p-3] + (*Problem).u_stag[i][(*Problem).Ny_p-4]);
    // Down
    for(int i = 1; i < (*Problem).Nx_p-2 ; i++ )
        (*Problem).u_stag[i][(*Problem).imap[i]+1] = -0.2*(15*(*Problem).u_stag[i][(*Problem).imap[i]+2] -5.0*(*Problem).u_stag[i][(*Problem).imap[i]+3] + (*Problem).u_stag[i][(*Problem).imap[i]+4]);
    // Side
    if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
        for(int j = 0; j < (*Problem).NHs-2; j++ )
            (*Problem).v_stag[(*Problem).NLs-1][j] =  -0.2*(15*(*Problem).v_stag[(*Problem).NLs-2][j] -5.0*(*Problem).v_stag[(*Problem).NLs-3][j] + (*Problem).v_stag[(*Problem).NLs-4][j]);

}

void boundary_psi_update(struct _problem* Problem, double (*Q)(struct _problem*)){
    // Upper boundary
    for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).psi[i][(*Problem).Ny-1] = Q(Problem);
    // Down boundary - Left - Right
    for(int i = 0; i < (*Problem).Nx ; i++ ){
        (*Problem).psi[i][(*Problem).imap[i]] = 0.0;
    }
    // Down boundary Side
    if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
        for(int j = 0; j < (*Problem).NHs; j++ ) (*Problem).psi[(*Problem).NLs-1][j] = 0.0;

}

void boundary_psi_save(struct _problem* Problem){
    for(int j = (*Problem).imap[0] ; j < (*Problem).Ny ; j++   )             (*Problem).psi_in[j] = (*Problem).psi[0][j];
    for(int j = (*Problem).imap[(*Problem).Nx-1] ; j < (*Problem).Ny ; j++ ) (*Problem).psi_out[j] = (*Problem).psi[0][j];
}

void boundary_omega_update(struct _problem* Problem){
    // Upper boundary
    for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][(*Problem).Ny-1   ] = -3.0/((*Problem).h*(*Problem).h) * ((*Problem).psi[i][(*Problem).Ny-2    ] - functionQ(Problem)) - 0.5*(*Problem).omega[i][(*Problem).Ny-2     ];
    // Down boundary - Left - Right
    for(int i = 0; i < (*Problem).Nx ; i++ ) (*Problem).omega[i][(*Problem).imap[i]] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[i][(*Problem).imap[i]+1] - 0.5*(*Problem).omega[i][(*Problem).imap[i]+1];
    // Down boundary Side
    if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
        for(int j = 0; j < (*Problem).NHs; j++ ) (*Problem).omega[(*Problem).NLs-1][j] = -3.0/((*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs+1][j] - 0.5*(*Problem).omega[(*Problem).NLs+1][j];
    // Corner boundary
    if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 ){
        (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][(*Problem).NHs] - 0.5* (*Problem).omega[(*Problem).NLs][(*Problem).NHs];
        (*Problem).omega[(*Problem).NLs-1][0] = - 3.0 / (2*(*Problem).h*(*Problem).h) * (*Problem).psi[(*Problem).NLs][1] - 0.5* (*Problem).omega[(*Problem).NLs][1];
    }

}

void boundary_omega_dwdx_update(struct _problem* Problem){
    // Inflow Boundary - Natural Condition
    for(int j = (*Problem).imap[0]; j < (*Problem).Ny-1 ; j++) (*Problem).omega[0][j] = (*Problem).omega[1][j];
    // Outflow Boundary - Natural Condition
    for(int i = (*Problem).imap[(*Problem).Nx-1]; i < (*Problem).Ny-1 ; i++) (*Problem).omega[(*Problem).Nx-1][i] = (*Problem).omega[(*Problem).Nx-2][i];
}
void boundary_pression_ghost_update(struct _problem* Problem){
    // Upper
    for(int i = 1; i < (*Problem).Nx_p-2 ; i++ )
        (*Problem).P[i][(*Problem).Ny_p-2] = (*Problem).P[i][(*Problem).Ny_p-3] + (*Problem).nu * ( (*Problem).omega[i][(*Problem).Ny-1] - (*Problem).omega[i-1][(*Problem).Ny-1] );
    // Down
    for(int i = 1; i < (*Problem).Nx_p-2 ; i++ )
        (*Problem).P[i][(*Problem).imap[i]] = (*Problem).P[i][(*Problem).imap[i]+1] - (*Problem).nu * ( (*Problem).omega[i][(*Problem).imap[i]] - (*Problem).omega[i-1][(*Problem).imap[i]+1] );
    // Side
    if( (*Problem).Ls != (*Problem).L && (*Problem).Ls != 0.0 )
        for(int j = 1; j < (*Problem).NHs; j++ )
            (*Problem).P[(*Problem).NLs-1][j] = (*Problem).P[(*Problem).NLs][j] + (*Problem).nu * ( (*Problem).omega[(*Problem).NLs-1][j] - (*Problem).omega[(*Problem).NLs-1][j-1] );
}
void boundary_pression_ghost_in_out(struct _problem* Problem){
    // Inflow

    for(int j = (*Problem).imap[0]+1; j < (*Problem).Ny_p-2 ; j++)
        (*Problem).P[0][j] =   (*Problem).P[1][j]
        + (((*Problem).psi[0][j]-(*Problem).psi[0][j-1]) - ((*Problem).psi_in[j]-(*Problem).psi_in[j-1]))/((*Problem).dt)
        + (*Problem).nu*((*Problem).omega[0][j]-(*Problem).omega[0][j-1]);
    // outflow

    for(int j = (*Problem).imap[(*Problem).Nx-1]+1; j < (*Problem).Ny_p-2 ; j++)
        (*Problem).P[(*Problem).Nx_p-2][j] = -(  (*Problem).P[(*Problem).Nx_p-3][j]
                                               + (( (*Problem).psi[(*Problem).Nx-1][j] - (*Problem).psi[(*Problem).Nx-1][j-1] ) - ((*Problem).psi_out[j]-(*Problem).psi_out[j-1]))/((*Problem).dt)
                                               + (*Problem).nu*((*Problem).omega[(*Problem).Nx-1][j]-(*Problem).omega[(*Problem).Nx-1][j-1]));

}
void boundary_pression_ghost_corner(struct _problem* Problem,int flag){
    // 1 - Right, else - Left
    // Vertical
    if(flag == 1) (*Problem).P[(*Problem).NLs-1][(*Problem).NHs-1] =
        (*Problem).P[(*Problem).NLs][(*Problem).NHs-1]  + (*Problem).nu * ( (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] - (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-2] );
    // Horizontal
    else          (*Problem).P[(*Problem).NLs-1][(*Problem).NHs-1] =
        (*Problem).P[(*Problem).NLs-1][(*Problem).NHs]  - (*Problem).nu * ( (*Problem).omega[(*Problem).NLs-1][(*Problem).NHs-1] - (*Problem).omega[(*Problem).NLs-2][(*Problem).NHs-1] );
}

void inner_psi_update(struct _problem* Problem){

    for(int i = 1; i < (*Problem).Nx-1; i++){
        for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
            (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
        }
    }/*
      for(int i = (*Problem).Nx-2; i > 0 ; i--){
      for(int j = 1; j < (*Problem).imax_map[i]-1; j++){
      (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
      }
    }*/ // Starting from another side - May be a potential improvement to switch side to ensure a faster convergence
}
void inner_pres_update(struct _problem* Problem){
    boundary_pression_ghost_corner(Problem,0);
    for(int i = 1; i < (*Problem).NLs; i++){
        for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny_p-2; j++){
            (*Problem).P[i][j] = scalar_pres_compute(Problem,i,j);
        }
    }
    boundary_pression_ghost_corner(Problem,1);
    for(int i =  (*Problem).NLs; i < (*Problem).Nx_p-2; i++){
        for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny_p-2; j++){
            (*Problem).P[i][j] = scalar_pres_compute(Problem,i,j);
        }
    }

    /*
     for(int i = (*Problem).Nx-2; i > 0 ; i--){
     for(int j = 1; j < (*Problem).imax_map[i]-1; j++){
     (*Problem).psi[i][j] = scalar_psi_compute(Problem,i,j);
     }
   }*/ // Starting from another side - May be a potential improvement to switch side to ensure a faster convergence
}

double inner_psi_error_compute(struct _problem* Problem){
    double e_error = 0.0;
    double square = 0.0;
    for(int i=1; i < (*Problem).Nx -1; i++){
        for(int j= (*Problem).imap[i]+1; j < (*Problem).Ny-1; j++){
            (*Problem).R_res[i][j] = scalar_psi_r_compute(Problem,i,j);
            square = (*Problem).R_res[i][j]*(*Problem).R_res[i][j] + square;
        }
    }
    e_error = fabs((*Problem).H*(*Problem).H/(*Problem).Q0*(*Problem).h*sqrt(1.0/((*Problem).L*(*Problem).H) *square));
    return e_error;
}
double inner_pres_error_compute(struct _problem* Problem){
    double square_L  = 0.0;
    double square_R  = 0.0;
    double e_error_L = 0.0;
    double e_error_R = 0.0;
    boundary_pression_ghost_corner(Problem,0);
    for(int i=1; i < (*Problem).NLs; i++){
        for(int j= (*Problem).imap[i]+1; j < (*Problem).Ny_p-2; j++){
            (*Problem).R_res_pres[i][j] = scalar_pres_r_compute(Problem,i,j);
            square_L = (*Problem).R_res_pres[i][j]*(*Problem).R_res_pres[i][j] + square_L;
        }
    }
    e_error_L = fabs((*Problem).H*(*Problem).H/(*Problem).Q0*(*Problem).h*sqrt(1.0/((*Problem).L*(*Problem).H) *square_L));
    boundary_pression_ghost_corner(Problem,1);
    for(int i =  (*Problem).NLs; i < (*Problem).Nx_p-2; i++){
        for(int j = (*Problem).imap[i]+1; j < (*Problem).Ny_p-2; j++){
            (*Problem).R_res_pres[i][j] = scalar_pres_r_compute(Problem,i,j);
            square_R = (*Problem).R_res_pres[i][j]*(*Problem).R_res_pres[i][j] + square_R;
        }
    }
    e_error_R = fabs((*Problem).H*(*Problem).H/(*Problem).Q0*(*Problem).h*sqrt(1.0/((*Problem).L*(*Problem).H) *square_R));
    return e_error_L+e_error_R;
}

void poisson_inner_psi_iterator(struct _problem* Problem){
    int n_iter = 0;
    int iter = 0 ;
    int iter_max = 5000;
    double error = (*Problem).tol+1;

    while( error>(*Problem).tol && iter < iter_max){
        n_iter++;
        inner_psi_update(Problem);
        // Inflow Boundary - Natural Condition
        for(int j = (*Problem).imap[0]; j < (*Problem).Ny-1 ; j++) (*Problem).psi[0][j]               = (*Problem).psi[1][j];
        // Outflow Boundary - Natural Condition
        for(int i = (*Problem).imap[(*Problem).Nx-1]; i < (*Problem).Ny -1 ; i++) (*Problem).psi[(*Problem).Nx-1][i] = (*Problem).psi[(*Problem).Nx-2][i];

        error = inner_psi_error_compute(Problem);
        iter++;
    }
    if(iter >= iter_max){fprintf(stderr, "[DEADSTOP] Maximum Psi Iteration Reached Current Error: %f\n",error); deadstop_exit(Problem);}
}
void poisson_inner_pres_iterator(struct _problem* Problem){
    int n_iter = 0;
    int iter = 0 ;
    int iter_max = 200000;
    double error = (*Problem).tol+1;

    while( error>(*Problem).tol && iter < iter_max){
        n_iter++;

        inner_pres_update(Problem);
        boundary_pression_ghost_in_out(Problem);
        boundary_pression_ghost_update(Problem);

        error = inner_pres_error_compute(Problem);
        iter++;
    }
    if(iter >= iter_max){fprintf(stderr, "[DEADSTOP] Maximum Pressure Iteration Reached Current Error: %f\n",error); deadstop_exit(Problem);}
}

int diagnose_check(struct _problem* Problem,int i,int j,int ktime){
    int check = 0;
    double Re_h       = (*Problem).h*(fabs((*Problem).u[i][j]) + fabs((*Problem).v[i][j]))/(*Problem).nu;
    if(Re_h > (*Problem).Re_h[ktime]) (*Problem).Re_h[ktime] = Re_h;
    double Re_h_omega = (*Problem).h*(*Problem).h*fabs((*Problem).omega[i][j])/(*Problem).nu;
    if(Re_h_omega > (*Problem).Re_h_omega[ktime]) (*Problem).Re_h_omega[ktime] = Re_h_omega;
    double Beta       = (fabs((*Problem).u[i][j]) + fabs((*Problem).v[i][j]))* (*Problem).dt/(*Problem).h;
    if(Beta > (*Problem).Beta_CFL[ktime]) (*Problem).Beta_CFL[ktime] = Beta;
    if( Re_h       >=  5 )                  check += 1;
    if( Re_h_omega >= 20 )                  check += 2;
    if( Beta       >=  (*Problem).CFL    )  check += 4;
    return check;
}
