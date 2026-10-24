/**********************************************************************
  Local_physical_quantities.c:

     Local_physical_quantities.c is a subrutine to output values of
     local physical quantities.

  Log of Local_physical_quantities.c:

     01/May/2015  Released by M.Fukuda
     13/Dec/2019  bug in kinetic energy
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"
#include "lapack_prototype.h"
#include "mpi.h"
#include <omp.h>
#include "Local_physical_quantities.h"

/*  scalar  */
static void electron_density(double *out, int BN);
static void charge_density(double *out, int BN);
static void zeta_potential(double *out, int BN);
static void kinetic_energy_density(double *out, int BN);
//static void energy_density(FILE *fp, int BN, int flag_cube);
static void charge_density_alpha(double *out, int BN);   /* added by yoshida 23.December.2019 */
static void charge_density_beta(double *out, int BN);   /* added by yoshida 23.December.2019 */
static void energy(double *out, int BN);   /* added by takamatsu 19.February.2020 */
static void chemical_potential(double *out, int BN);   /* added by takamatsu 19.February.2020 */

/*  vector  */
static void charge_current(double *out, int BN);
static void kinetic_momentum(double *out, int BN);
static void kinetic_momentum_alpha(double *out, int BN);
static void kinetic_momentum_beta(double *out, int BN);
//static void kinetic_momentum_spinup(FILE *fp, int BN);
//static void kinetic_momentum_spindown(FILE *fp, int BN);
static void spin_vorticity(double *out, int BN);
static void spin_angular_momentum(double *out, int BN);
//static void spin_torque(FILE *fp, int BN);
static void zeta_force(double *out, int BN);
static void pau_dd(double *out, int BN);
//static void tension(FILE *fp, int BN);
//static void lorentz_force(FILE *fp, int BN);
static void gradient_electron_density(double *out, int BN);

static void spin_z_current(double *out, int BN);   /* added by kuroda 7.June.2019 */
static void charge_current_alpha(double *out, int BN);     /* added by kuroda 5.December.2019 */
static void charge_current_beta(double *out, int BN);   /* added by kuroda 5.December.2019 */
static void velocity(double *out, int BN);   /* added by kuroda 31.January.2020 */
static void torq_zeta(double *out, int BN);      /* added by shimada 27.December.2019 */
static void soc_part1(double *out, int BN);       /* added by shimada September.2020 */
static void soc_part2(double *out, int BN);       /* added by shimada September.2020 */

/*  tensor  */
static void stress(double *out, int BN);
static void compute_stress_full(double val[3][3], int BN);
static void stress_S(double *out, int BN);
static void stress_A(double *out, int BN);
static void stress_diag(double *out, int BN);
static void dd_electron_density(double *out, int BN);
static void stress_nonrel(double *out, int BN);   /* added by yoshida 7.January.2020 */
static void stress_nonrel_alpha(double *out, int BN);   /* added by yoshida 7.January.2020 */
static void stress_nonrel_beta(double *out, int BN);   /* added by yoshida 7.January.2020 */

/* diagonalization */
static void Eigen_lapack_d(double **a, double *ko, int n, int EVmax);

/* In-memory values use row-major Cartesian tensor ordering. */
static void store_scalar(double *out, double val) { out[0] = val; }
static void store_vector(double *out, const double val[3])
{
  int i;
  for (i=0; i<3; i++) out[i] = val[i];
}
static void store_tensor(double *out, double val[3][3])
{
  int i,j;
  for (i=0; i<3; i++) for (j=0; j<3; j++) out[3*i+j] = val[i][j];
}
static void store_tensor_diag(double *out, double val[3][3], const double ko[3])
{
  int i;
  for (i=0; i<3; i++) out[i] = ko[i];
  store_tensor(out+3,val);
}

/* print */
static void print_scalar(FILE *fp, double val, int flag_cube);
static void print_vector(FILE *fp, double val[3]);
static void print_tensor(FILE *fp, double val[3][3]);
static void print_tensor_diag(FILE *fp, double val[3][3], double ko[3]);

//char *FMT1 = "%24.14E";
char *FMT1 = "%18.8E";
char *FMT_cube = "%13.3E";

/******************************************
                print format
 ******************************************/

void print_scalar(FILE *fp, double val, int flag_cube)
{
    if(flag_cube==1){
      fprintf(fp,FMT_cube,val);
    }
    else {
      fprintf(fp,FMT1,val);
      fprintf(fp,"\n");
    }
}

void print_vector(FILE *fp, double val[3])
{
   int i;
   for (i=0; i<3; i++){
      fprintf(fp,FMT1,val[i]);
   }
   fprintf(fp,"\n");
}

void print_tensor(FILE *fp, double val[3][3])
{
   int i,j;
   for (j=0; j<3; j++){
      for (i=0; i<3; i++){
       fprintf(fp,FMT1, val[i][j]);
      }
   }
   fprintf(fp,"\n");
}

void print_tensor_diag(FILE *fp, double val[3][3], double ko[3])
{
   int i,j;
   for (i=0; i<3; i++){
      fprintf(fp,FMT1, ko[i]);
   }
   for (j=0; j<3; j++){
      for (i=0; i<3; i++){
       fprintf(fp,FMT1, val[i][j]);
      }
   }
   fprintf(fp,"\n");
}




/******************************************
             selection part
 ******************************************/

/* Evaluate local physical quantities independently of output formatting. */
int flpq_evaluate_local_quantity(int rank, int n, int BN,
                                  FLPQ_LocalValue *result)
{
  if (result == NULL) return 0;
  memset(result, 0, sizeof(*result));
  result->rank = rank;
  result->index = n;
  switch (rank) {
    case 0:
      result->components = 1;
      switch (n) {
        case 0: electron_density(result->values, BN); return 1;
        case 1: charge_density(result->values, BN); return 1;
        case 2: zeta_potential(result->values, BN); return 1;
        case 3: kinetic_energy_density(result->values, BN); return 1;
        case 5: charge_density_alpha(result->values, BN); return 1;
        case 6: charge_density_beta(result->values, BN); return 1;
        case 7: energy(result->values, BN); return 1;
        case 8: chemical_potential(result->values, BN); return 1;
        default: return 0;
      }
    case 1:
      result->components = 3;
      switch (n) {
        case 0: charge_current(result->values, BN); return 1;
        case 1: kinetic_momentum(result->values, BN); return 1;
        case 2: spin_vorticity(result->values, BN); return 1;
        case 3: spin_angular_momentum(result->values, BN); return 1;
        case 5: zeta_force(result->values, BN); return 1;
        case 6: pau_dd(result->values, BN); return 1;
        case 9: gradient_electron_density(result->values, BN); return 1;
        case 10: spin_z_current(result->values, BN); return 1;
        case 11: charge_current_alpha(result->values, BN); return 1;
        case 12: charge_current_beta(result->values, BN); return 1;
        case 13: velocity(result->values, BN); return 1;
        case 14: torq_zeta(result->values, BN); return 1;
        case 15: soc_part1(result->values, BN); return 1;
        case 16: kinetic_momentum_alpha(result->values, BN); return 1;
        case 17: kinetic_momentum_beta(result->values, BN); return 1;
        case 18: soc_part2(result->values, BN); return 1;
        default: return 0;
      }
    case 2:
      result->components = (n == 3 ? 12 : 9);
      switch (n) {
        case 0: stress(result->values, BN); return 1;
        case 1: stress_S(result->values, BN); return 1;
        case 2: stress_A(result->values, BN); return 1;
        case 3: stress_diag(result->values, BN); return 1;
        case 4: dd_electron_density(result->values, BN); return 1;
        case 5: stress_nonrel(result->values, BN); return 1;
        case 6: stress_nonrel_alpha(result->values, BN); return 1;
        case 7: stress_nonrel_beta(result->values, BN); return 1;
        default: return 0;
      }
    default: return 0;
  }
}

/* Preserve legacy numeric formats and tensor column ordering. */
void flpq_write_local_value(FILE *fp, const FLPQ_LocalValue *result,
                            int flag_cube)
{
  double matrix[3][3];
  int i, j;
  if (fp == NULL || result == NULL) return;
  if (result->rank == 0) print_scalar(fp, result->values[0], flag_cube);
  else if (result->rank == 1) print_vector(fp, (double *)result->values);
  else if (result->rank == 2) {
    int off = (result->index == 3) ? 3 : 0;
    for (i=0; i<3; i++) for (j=0; j<3; j++)
      matrix[i][j] = result->values[off+3*i+j];
    if (result->index == 3) {
      double eigenvalues[3];
      for (i=0; i<3; i++) eigenvalues[i] = result->values[i];
      print_tensor_diag(fp, matrix, eigenvalues);
    }
    else print_tensor(fp, matrix);
  }
}

void scalar_density(int n, FILE *fp, int BN, int flag_cube)
{
  FLPQ_LocalValue v;
  if (flpq_evaluate_local_quantity(0,n,BN,&v))
    flpq_write_local_value(fp,&v,flag_cube);
}
void vector_density(int n, FILE *fp, int BN)
{
  FLPQ_LocalValue v;
  if (flpq_evaluate_local_quantity(1,n,BN,&v))
    flpq_write_local_value(fp,&v,0);
}
void tensor_density(int n, FILE *fp, int BN)
{
  FLPQ_LocalValue v;
  if (flpq_evaluate_local_quantity(2,n,BN,&v))
    flpq_write_local_value(fp,&v,0);
}

/******************************************
        Local physical quantities
 ******************************************/
/*******************************
             scalar
 *******************************/
void electron_density(double *out, int BN)
{ 
    double val;

    if (SpinP_switch==0) {
       val = 2.0*Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN] + Density_Grid_B[1][BN];
    }
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);

}

void charge_density(double *out, int BN)
{ 
    double val;

    if (SpinP_switch==0) {
       val = 2.0*Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN] + Density_Grid_B[1][BN];
    }
    val = - val;
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);

}

void zeta_potential(double *out, int BN)
{ 
    double val;

    val =  dDensity_Grid_B[0][3][BN] - dDensity_Grid_B[0][5][BN]
          -dDensity_Grid_B[1][2][BN] + dDensity_Grid_B[1][4][BN]
          +dDensity_Grid_B_i[2][0][BN] - dDensity_Grid_B_i[2][1][BN];
    val = 0.5*val;
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);

}

void kinetic_energy_density(double *out, int BN)
{ 
    double val;

    //val =  ddDensity_Grid_B[0][0][BN] + ddDensity_Grid_B[0][1][BN];
    //val = -0.5*val;
    if (SpinP_switch==0) {
       val =       ddDensity_Grid_B[1][0][BN];
       val = val + ddDensity_Grid_B[2][0][BN];
       val = val + ddDensity_Grid_B[3][0][BN];
       val = -val;
    }
    else {
       val =       ddDensity_Grid_B[1][0][BN] + ddDensity_Grid_B[1][1][BN];
       val = val + ddDensity_Grid_B[2][0][BN] + ddDensity_Grid_B[2][1][BN];
       val = val + ddDensity_Grid_B[3][0][BN] + ddDensity_Grid_B[3][1][BN];
       val = -0.5*val;
    }
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);

}


void charge_density_alpha(double *out, int BN)   /* added by yoshida 23.December.2019 */
{ 
    double val;

    if (SpinP_switch==0) {
       val = Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN];
    }
    val = - val;
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);
}

void charge_density_beta(double *out, int BN)   /* added by yoshida 23.December.2019 */
{ 
    double val;

    if (SpinP_switch==0) {
       val = Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[1][BN];
    }
    val = - val;
    store_scalar(out, val);
    //printf("%24.14E %d\n", val, BN);
}

void energy(double *out, int BN)
{
    double val[3][3];
    compute_stress_full(val,BN);
    /* Eigenvalue sum of the symmetric stress equals its trace. */
    store_scalar(out, 0.5*(val[0][0]+val[1][1]+val[2][2]));
}

void chemical_potential(double *out, int BN)
{
    double val[3][3], density, energy_density, value;
    compute_stress_full(val,BN);
    energy_density=0.5*(val[0][0]+val[1][1]+val[2][2]);
    if (SpinP_switch==0)
      density=2.0*Density_Grid_B[0][BN];
    else
      density=Density_Grid_B[0][BN]+Density_Grid_B[1][BN];
    /* Retain the original regularization at vanishing density. */
    value=fabs(density)<1.e-16 ? 0.0 : energy_density/density;
    store_scalar(out, value);
}

/*******************************
            vector 
 *******************************/
void charge_current(double *out, int BN)
{
  /* minus signs was modified in 3.Aug.2015 */
    double val[3];
    
    if (SpinP_switch==0) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = 0.0;

       if (Solver==4) {
          val[0] = -2.0*dDensity_Grid_B_i[0][0][BN];
          val[1] = -2.0*dDensity_Grid_B_i[1][0][BN];
          val[2] = -2.0*dDensity_Grid_B_i[2][0][BN];
       }
    }
    else if (SpinP_switch==1) {
       val[0] = 0.5*(-dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
       val[1] = 0.5*( dDensity_Grid_B[0][0][BN] - dDensity_Grid_B[0][1][BN] );
       val[2] = 0.0;
       if (Solver==4) {
         val[0] = val[0] -( dDensity_Grid_B_i[0][0][BN] + dDensity_Grid_B_i[0][1][BN] );
         val[1] = val[1] -( dDensity_Grid_B_i[1][0][BN] + dDensity_Grid_B_i[1][1][BN] );
         val[2] = val[2] -( dDensity_Grid_B_i[2][0][BN] + dDensity_Grid_B_i[2][1][BN] );
       }
    }
    else if (SpinP_switch==3) {
       val[0] = 0.5*( dDensity_Grid_B[2][3][BN] + dDensity_Grid_B[2][5][BN] 
                     -dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
       val[1] = 0.5*(-dDensity_Grid_B[2][2][BN] - dDensity_Grid_B[2][4][BN] 
                     +dDensity_Grid_B[0][0][BN] - dDensity_Grid_B[0][1][BN] );
       val[2] = 0.5*( dDensity_Grid_B[1][2][BN] + dDensity_Grid_B[1][4][BN] 
                     -dDensity_Grid_B[0][3][BN] - dDensity_Grid_B[0][5][BN] );
    
       val[0] = val[0] -( dDensity_Grid_B_i[0][0][BN] + dDensity_Grid_B_i[0][1][BN] );
       val[1] = val[1] -( dDensity_Grid_B_i[1][0][BN] + dDensity_Grid_B_i[1][1][BN] );
       val[2] = val[2] -( dDensity_Grid_B_i[2][0][BN] + dDensity_Grid_B_i[2][1][BN] );
    }
    store_vector(out, val);
}

void kinetic_momentum(double *out, int BN)
{
    double val[3];
    
    if (SpinP_switch==0) {
      val[0] = 2.0*dDensity_Grid_B_i[0][0][BN];
      val[1] = 2.0*dDensity_Grid_B_i[1][0][BN];
      val[2] = 2.0*dDensity_Grid_B_i[2][0][BN];
    }
    else {
      val[0] = dDensity_Grid_B_i[0][0][BN] + dDensity_Grid_B_i[0][1][BN];
      val[1] = dDensity_Grid_B_i[1][0][BN] + dDensity_Grid_B_i[1][1][BN];
      val[2] = dDensity_Grid_B_i[2][0][BN] + dDensity_Grid_B_i[2][1][BN];
    }
    store_vector(out, val);
}

void kinetic_momentum_alpha(double *out, int BN)
{
    double val[3];
    
    val[0] = dDensity_Grid_B_i[0][0][BN];
    val[1] = dDensity_Grid_B_i[1][0][BN];
    val[2] = dDensity_Grid_B_i[2][0][BN];
    store_vector(out, val);
}

void kinetic_momentum_beta(double *out, int BN)
{
    double val[3];
    
    val[0] = dDensity_Grid_B_i[0][1][BN];
    val[1] = dDensity_Grid_B_i[1][1][BN];
    val[2] = dDensity_Grid_B_i[2][1][BN];
    store_vector(out, val);
}

//void kinetic_momentum_spinup(FILE *fp, int BN)
//{
//    double val[3];
//    
//    val[0] = dDensity_Grid_B_i[0][0][BN];
//    val[1] = dDensity_Grid_B_i[1][0][BN];
//    val[2] = dDensity_Grid_B_i[2][0][BN];
//
//    print_vector(fp, val);
//}
//
//void kinetic_momentum_spindown(FILE *fp, int BN)
//{
//    double val[3];
//    
//    val[0] = dDensity_Grid_B_i[0][1][BN];
//    val[1] = dDensity_Grid_B_i[1][1][BN];
//    val[2] = dDensity_Grid_B_i[2][1][BN];
//
//    print_vector(fp, val);
//}

void spin_vorticity(double *out, int BN)
{
  /* minus signs was modified in 3.Aug.2015 */
    double val[3];
    
    if (SpinP_switch==0) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = 0.0;
    }
    else if (SpinP_switch==1) {
       val[0] = -0.5*(-dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
       val[1] = -0.5*( dDensity_Grid_B[0][0][BN] - dDensity_Grid_B[0][1][BN] );
       val[2] =  0.0;
    }
    else if ((SpinP_switch==3)||(Solver==4)) {
       val[0] = -0.5*( dDensity_Grid_B[2][3][BN] + dDensity_Grid_B[2][5][BN] 
                      -dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
       val[1] = -0.5*(-dDensity_Grid_B[2][2][BN] - dDensity_Grid_B[2][4][BN] 
                      +dDensity_Grid_B[0][0][BN] - dDensity_Grid_B[0][1][BN] );
       val[2] = -0.5*( dDensity_Grid_B[1][2][BN] + dDensity_Grid_B[1][4][BN] 
                      -dDensity_Grid_B[0][3][BN] - dDensity_Grid_B[0][5][BN] );
    }
    store_vector(out, val);
}

void spin_angular_momentum(double *out, int BN)
{
    double val[3];
    
    if (SpinP_switch==0) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = 0.0;
    }
    else if (SpinP_switch==1) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = (Density_Grid_B[0][BN] - Density_Grid_B[1][BN])*0.5 ;
    }
    else if ((SpinP_switch==3)||(Solver==4)) {
       val[0] = Density_Grid_B[2][BN];
       val[1] = Density_Grid_B[3][BN];
       val[2] = (Density_Grid_B[0][BN] - Density_Grid_B[1][BN])*0.5 ;
    }
    store_vector(out, val);
}

//void spin_torque(FILE *fp, int BN)
//{
//    double val[3];
//    
//    val[0] = 0.0;
//    val[1] = 0.0;
//    val[2] = 0.0;
//    print_vector(fp, val);
//}
//
//

void torq_zeta(double *out, int BN) /* added by matsumoto shimada 27.December.2019 */
{ 
  int i,j;
  double a1[3], a2[3];
  double val[3][3];

  // double t0=0, t1=0, t2=0, z0=0, z1=0, z2=0;
  // spin_torque(fp, BN, t0, t1, t2);
  // zeta_force(fp, BN, z0, z1, z2);

  // val[0]=t0+z0; 
  // val[1]=t1+z1;
  // val[2]=t2+z2;

  // spin torque
//start of calculation of strss_A----------------------------------------------------------------------
    if (SpinP_switch==0) {
       val[0][0] = 0.0;
       val[1][1] = 0.0;
       val[2][2] = 0.0;
       val[0][1] = -( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==1) {
       val[0][0] = 0.0; 
       val[1][1] = 0.0; 
       val[2][2] = 0.0; 
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                        + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                        + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                        + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==3) {
       val[0][0] = 0.0; 
       val[1][1] = 0.0; 
       val[2][2] = 0.0; 
       val[0][1] =   dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                   + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] ;
       val[1][2] =   dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                   + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] ;
       val[2][0] =   dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                   + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] ;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

       val[0][1] = val[0][1] - dd2Density_Grid_B[6][5][BN] - ddDensity_Grid_B[6][5][BN]
                             + dd2Density_Grid_B[6][3][BN] + ddDensity_Grid_B[6][3][BN];
       val[1][0] = val[1][0] - dd2Density_Grid_B[5][2][BN] - ddDensity_Grid_B[5][4][BN]
                             + dd2Density_Grid_B[5][4][BN] + ddDensity_Grid_B[5][2][BN];
       val[1][2] = val[1][2] + dd2Density_Grid_B[2][5][BN] + ddDensity_Grid_B[2][5][BN]
                             - dd2Density_Grid_B[2][3][BN] - ddDensity_Grid_B[2][3][BN]
                             + dd2Density_Grid_B[4][4][BN] + ddDensity_Grid_B[4][4][BN]
                             - dd2Density_Grid_B[4][2][BN] - ddDensity_Grid_B[4][2][BN];
       val[2][1] = val[2][1] - dd2Density_Grid_B[3][5][BN] - ddDensity_Grid_B[3][5][BN]
                             + dd2Density_Grid_B[3][3][BN] + ddDensity_Grid_B[3][3][BN];
       val[2][0] = val[2][0] - dd2Density_Grid_B[3][4][BN] - ddDensity_Grid_B[3][4][BN]
                             + dd2Density_Grid_B[3][2][BN] + ddDensity_Grid_B[3][2][BN];
       val[0][2] = val[0][2] + dd2Density_Grid_B[4][3][BN] + ddDensity_Grid_B[4][5][BN]
                             - dd2Density_Grid_B[4][5][BN] - ddDensity_Grid_B[4][3][BN]
                             + dd2Density_Grid_B[1][4][BN] + ddDensity_Grid_B[1][4][BN]
                             - dd2Density_Grid_B[1][2][BN] - ddDensity_Grid_B[1][2][BN];
       /* for NEGF */
       val[0][0] = val[0][0] - dd2Density_Grid_B_i[4][0][BN] + ddDensity_Grid_B_i[4][0][BN]
                             + dd2Density_Grid_B_i[4][1][BN] - ddDensity_Grid_B_i[4][1][BN];
       val[1][1] = val[1][1] - dd2Density_Grid_B_i[4][0][BN] - ddDensity_Grid_B_i[4][0][BN]
                             + dd2Density_Grid_B_i[4][1][BN] + ddDensity_Grid_B_i[4][1][BN];
       val[0][1] = val[0][1] - dd2Density_Grid_B_i[1][0][BN] - ddDensity_Grid_B_i[1][0][BN]
                             + dd2Density_Grid_B_i[1][1][BN] + ddDensity_Grid_B_i[1][1][BN];
       val[1][0] = val[1][0] + dd2Density_Grid_B_i[2][0][BN] + ddDensity_Grid_B_i[2][0][BN]
                             - dd2Density_Grid_B_i[2][1][BN] - ddDensity_Grid_B_i[2][1][BN];
       val[2][1] = val[2][1] + dd2Density_Grid_B_i[6][0][BN] - ddDensity_Grid_B_i[6][0][BN]
                             - dd2Density_Grid_B_i[6][1][BN] + ddDensity_Grid_B_i[6][1][BN];
       val[2][0] = val[2][0] + dd2Density_Grid_B_i[5][0][BN] + ddDensity_Grid_B_i[5][0][BN]
                             - dd2Density_Grid_B_i[5][1][BN] - ddDensity_Grid_B_i[5][1][BN];

       /* anti-symmetrize */
       val[0][1] = (val[0][1] - val[1][0])*0.5;
       val[1][2] = (val[1][2] - val[2][1])*0.5;
       val[2][0] = (val[2][0] - val[0][2])*0.5;
       val[1][0] = - val[0][1];
       val[2][1] = - val[1][2];
       val[0][2] = - val[2][0];

       val[0][0] = 0.0;   /* added by yoshida 25.December.2019 */
       val[1][1] = 0.0;   /* added by yoshida 25.December.2019 */ 
       val[2][2] = 0.0;   /* added by yoshida 25.December.2019 */

       /* coef -0.5 */
       for(i=0;i<3;i++){
           for(j=0;j<3;j++){
               val[i][j] = -0.5*val[i][j];
           }
       }

    }
//end of calculation of strss_A----------------------------------------------------------------------
   
    a1[0] = 2.0*val[2][1] ; /* added by shimada 27.December.2019 */
    a1[1] = 2.0*val[0][2] ; /* added by shimada 27.December.2019 */ 
    a1[2] = 2.0*val[1][0] ; /* added by shimada 27.December.2019 */

    // zeta force
    a2[0]  = -0.5*( ddDensity_Grid_B[1][3][BN] + dd2Density_Grid_B[1][3][BN] 
                   -ddDensity_Grid_B[1][5][BN] - dd2Density_Grid_B[1][5][BN]
                   -ddDensity_Grid_B[4][2][BN] - dd2Density_Grid_B[4][2][BN] 
                   +ddDensity_Grid_B[4][4][BN] + dd2Density_Grid_B[4][4][BN]
                   +ddDensity_Grid_B_i[6][0][BN] - dd2Density_Grid_B_i[6][0][BN] 
                   -ddDensity_Grid_B_i[6][1][BN] + dd2Density_Grid_B_i[6][1][BN] );

    a2[1]  = -0.5*( ddDensity_Grid_B[4][3][BN] + dd2Density_Grid_B[4][5][BN] 
                   -ddDensity_Grid_B[4][5][BN] - dd2Density_Grid_B[4][3][BN]
                   -ddDensity_Grid_B[2][2][BN] - dd2Density_Grid_B[2][2][BN] 
                   +ddDensity_Grid_B[2][4][BN] + dd2Density_Grid_B[2][4][BN]
                   +ddDensity_Grid_B_i[5][0][BN] + dd2Density_Grid_B_i[5][0][BN] 
                   -ddDensity_Grid_B_i[5][1][BN] - dd2Density_Grid_B_i[5][1][BN] );

    a2[2]  = -0.5*( ddDensity_Grid_B[6][3][BN] + dd2Density_Grid_B[6][3][BN] 
                   -ddDensity_Grid_B[6][5][BN] - dd2Density_Grid_B[6][5][BN]
                   -ddDensity_Grid_B[5][2][BN] - dd2Density_Grid_B[5][4][BN] 
                   +ddDensity_Grid_B[5][4][BN] + dd2Density_Grid_B[5][2][BN]
                   +ddDensity_Grid_B_i[3][0][BN] + dd2Density_Grid_B_i[3][0][BN] 
                   -ddDensity_Grid_B_i[3][1][BN] - dd2Density_Grid_B_i[3][1][BN] );

    // spin torque + zeta force

    a1[0] = a1[0] + a2[0];
    a1[1] = a1[1] + a2[1];
    a1[2] = a1[2] + a2[2];
    store_vector(out, a1);
}
void zeta_force(double *out, int BN)
{
    double val[3];
    
    val[0] = -0.5*( ddDensity_Grid_B[1][3][BN] + dd2Density_Grid_B[1][3][BN] 
                   -ddDensity_Grid_B[1][5][BN] - dd2Density_Grid_B[1][5][BN]
                   -ddDensity_Grid_B[4][2][BN] - dd2Density_Grid_B[4][2][BN] 
                   +ddDensity_Grid_B[4][4][BN] + dd2Density_Grid_B[4][4][BN]
                   +ddDensity_Grid_B_i[6][0][BN] - dd2Density_Grid_B_i[6][0][BN] 
                   -ddDensity_Grid_B_i[6][1][BN] + dd2Density_Grid_B_i[6][1][BN] );

    val[1] = -0.5*( ddDensity_Grid_B[4][3][BN] + dd2Density_Grid_B[4][5][BN] 
                   -ddDensity_Grid_B[4][5][BN] - dd2Density_Grid_B[4][3][BN]
                   -ddDensity_Grid_B[2][2][BN] - dd2Density_Grid_B[2][2][BN] 
                   +ddDensity_Grid_B[2][4][BN] + dd2Density_Grid_B[2][4][BN]
                   +ddDensity_Grid_B_i[5][0][BN] + dd2Density_Grid_B_i[5][0][BN] 
                   -ddDensity_Grid_B_i[5][1][BN] - dd2Density_Grid_B_i[5][1][BN] );

    val[2] = -0.5*( ddDensity_Grid_B[6][3][BN] + dd2Density_Grid_B[6][3][BN] 
                   -ddDensity_Grid_B[6][5][BN] - dd2Density_Grid_B[6][5][BN]
                   -ddDensity_Grid_B[5][2][BN] - dd2Density_Grid_B[5][4][BN] 
                   +ddDensity_Grid_B[5][4][BN] + dd2Density_Grid_B[5][2][BN]
                   +ddDensity_Grid_B_i[3][0][BN] + dd2Density_Grid_B_i[3][0][BN] 
                   -ddDensity_Grid_B_i[3][1][BN] - dd2Density_Grid_B_i[3][1][BN] );
    store_vector(out, val);
}

void pau_dd(double *out, int BN)
{
    double val[3];

    val[0] = -0.5*( ddDensity_Grid_B[0][3][BN] - ddDensity_Grid_B[0][5][BN] );
    val[1] = -0.5*(-ddDensity_Grid_B[0][2][BN] + ddDensity_Grid_B[0][4][BN] );
    val[2] = -0.5*(-ddDensity_Grid_B_i[0][0][BN] - ddDensity_Grid_B_i[0][1][BN] );
    store_vector(out, val);
}

void gradient_electron_density(double *out, int BN)
{
    double val[3];

    val[0] = 2.0*( dDensity_Grid_B[0][0][BN] + dDensity_Grid_B[0][1][BN] );
    val[1] = 2.0*( dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
    val[2] = 2.0*( dDensity_Grid_B[2][0][BN] + dDensity_Grid_B[2][1][BN] );
    store_vector(out, val);
}


void spin_z_current(double *out, int BN)  /* added by kuroda 7.June.2019 */
{
    double val[3];
    
    val[0] = 0.5*(dDensity_Grid_B_i[0][0][BN] - dDensity_Grid_B_i[0][1][BN]);
    val[1] = 0.5*(dDensity_Grid_B_i[1][0][BN] - dDensity_Grid_B_i[1][1][BN]);
    val[2] = 0.5*(dDensity_Grid_B_i[2][0][BN] - dDensity_Grid_B_i[2][1][BN]);
    store_vector(out, val);
}


void charge_current_alpha(double *out, int BN)  /* added by kuroda 5.December.2019 */
{
    double val[3];
    
    val[0] = -dDensity_Grid_B_i[0][0][BN];
    val[1] = -dDensity_Grid_B_i[1][0][BN];
    val[2] = -dDensity_Grid_B_i[2][0][BN];
    store_vector(out, val);
}

void charge_current_beta(double *out, int BN)  /* added by kuroda 5.December.2019 */
{
    double val[3];
    
    val[0] = -dDensity_Grid_B_i[0][1][BN];
    val[1] = -dDensity_Grid_B_i[1][1][BN];
    val[2] = -dDensity_Grid_B_i[2][1][BN];
    store_vector(out, val);
}


void velocity(double *out, int BN)  /* added by kuroda 31.January.2020 */
{
    double val[3];
    
    val[0] =   dDensity_Grid_B_i[0][0][BN] + dDensity_Grid_B_i[0][1][BN] 
             + dDensity_Grid_B[1][0][BN] - dDensity_Grid_B[1][1][BN]
             - dDensity_Grid_B[2][5][BN] - dDensity_Grid_B[2][3][BN];

    val[1] =   dDensity_Grid_B_i[1][0][BN] + dDensity_Grid_B_i[1][1][BN] 
             + dDensity_Grid_B[2][4][BN] + dDensity_Grid_B[2][2][BN]
             - dDensity_Grid_B[0][0][BN] + dDensity_Grid_B[0][1][BN];

    val[2] =   dDensity_Grid_B_i[2][0][BN] + dDensity_Grid_B_i[2][1][BN] 
             + dDensity_Grid_B[0][5][BN] + dDensity_Grid_B[0][3][BN]
             - dDensity_Grid_B[0][4][BN] - dDensity_Grid_B[1][2][BN];

    store_vector(out, val);
}

void soc_part1(double *out, int BN) /* added by shimada September.2020 */
{
    double val[3];
    
    if (SpinP_switch==0) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = 0.0;
    }
    else if (SpinP_switch==1) {
       val[0] = (Density_Grid_B[2][BN])*0.25*0.00729735257*0.00729735257 ;
       val[1] = (Density_Grid_B[3][BN])*0.25*0.00729735257*0.00729735257 ;
       val[2] = (Density_Grid_B[0][BN] - Density_Grid_B[1][BN])*0.125*0.00729735257*0.00729735257 ;
    }
    else if (SpinP_switch==3) {

       val[0] = (Density_Grid_B[2][BN])*0.25*0.00729735257*0.00729735257 ;
       val[1] = (Density_Grid_B[3][BN])*0.25*0.00729735257*0.00729735257 ;
       val[2] = (Density_Grid_B[0][BN] - Density_Grid_B[1][BN])*0.125*0.00729735257*0.00729735257 ;
    }
    store_vector(out, val);
}

void soc_part2(double *out, int BN) /* added by shimada September.2020 */
{
    double val[3];
    
    if (SpinP_switch==0) {
       val[0] = 0.0;
       val[1] = 0.0;
       val[2] = 0.0;
    }
    else if (SpinP_switch==1) {
       val[0] = (Density_Grid_B[0][BN])*0.125*0.00729735257*0.00729735257 ;
       val[1] = (-Density_Grid_B[1][BN])*0.125*0.00729735257*0.00729735257 ;
       val[2] = 0.0 ;
    }
    else if (SpinP_switch==3) {

       val[0] = (Density_Grid_B[0][BN])*0.125*0.00729735257*0.00729735257 ;
       val[1] = (-Density_Grid_B[1][BN])*0.125*0.00729735257*0.00729735257 ;
       val[2] = 0.0 ;
    }
    store_vector(out, val);
}
/*******************************
             tensor
 *******************************/
static void compute_stress_full(double val[3][3], int BN)
{ 
    int i,j;

    if (SpinP_switch==0) {
       val[0][0] = -( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN] );
       val[1][1] = -( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN] );
       val[2][2] = -( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN] );
       val[0][1] = -( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==1) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN]
                        + dd2Density_Grid_B[1][1][BN] - ddDensity_Grid_B[1][1][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN]
                        + dd2Density_Grid_B[2][1][BN] - ddDensity_Grid_B[2][1][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN]
                        + dd2Density_Grid_B[3][1][BN] - ddDensity_Grid_B[3][1][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                        + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                        + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                        + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] );
        /* Initialize transposed entries before NEGF corrections. */
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
       /* for NEGF */
       if (Solver==4) {
         val[0][0] = val[0][0] - dd2Density_Grid_B_i[4][0][BN] + ddDensity_Grid_B_i[4][0][BN]
                               + dd2Density_Grid_B_i[4][1][BN] - ddDensity_Grid_B_i[4][1][BN];
         val[1][1] = val[1][1] - dd2Density_Grid_B_i[4][0][BN] - ddDensity_Grid_B_i[4][0][BN]
                               + dd2Density_Grid_B_i[4][1][BN] + ddDensity_Grid_B_i[4][1][BN];
         val[0][1] = val[0][1] - dd2Density_Grid_B_i[1][0][BN] - ddDensity_Grid_B_i[1][0][BN]
                               + dd2Density_Grid_B_i[1][1][BN] + ddDensity_Grid_B_i[1][1][BN];
         val[1][0] = val[1][0] + dd2Density_Grid_B_i[2][0][BN] + ddDensity_Grid_B_i[2][0][BN]
                               - dd2Density_Grid_B_i[2][1][BN] - ddDensity_Grid_B_i[2][1][BN];
         val[2][1] = val[2][1] + dd2Density_Grid_B_i[6][0][BN] - ddDensity_Grid_B_i[6][0][BN]
                               - dd2Density_Grid_B_i[6][1][BN] + ddDensity_Grid_B_i[6][1][BN];
         val[2][0] = val[2][0] + dd2Density_Grid_B_i[5][0][BN] + ddDensity_Grid_B_i[5][0][BN]
                               - dd2Density_Grid_B_i[5][1][BN] - ddDensity_Grid_B_i[5][1][BN];
       }

     }
    else if (SpinP_switch==3) {
       val[0][0] =   dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN]
                   + dd2Density_Grid_B[1][1][BN] - ddDensity_Grid_B[1][1][BN] ;
       val[1][1] =   dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN]
                   + dd2Density_Grid_B[2][1][BN] - ddDensity_Grid_B[2][1][BN] ;
       val[2][2] =   dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN]
                   + dd2Density_Grid_B[3][1][BN] - ddDensity_Grid_B[3][1][BN] ;
       val[0][1] =   dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                   + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] ;
       val[1][2] =   dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                   + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] ;
       val[2][0] =   dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                   + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] ;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

       val[0][0] = val[0][0] - dd2Density_Grid_B[6][4][BN] - ddDensity_Grid_B[6][4][BN]
                             + dd2Density_Grid_B[6][2][BN] + ddDensity_Grid_B[6][2][BN];
       val[1][1] = val[1][1] - dd2Density_Grid_B[5][3][BN] - ddDensity_Grid_B[5][5][BN] 
                             + dd2Density_Grid_B[5][5][BN] + ddDensity_Grid_B[5][3][BN];
       val[2][2] = val[2][2] + dd2Density_Grid_B[5][5][BN] + ddDensity_Grid_B[5][5][BN]
                             - dd2Density_Grid_B[5][3][BN] - ddDensity_Grid_B[5][3][BN]
                             + dd2Density_Grid_B[6][2][BN] + ddDensity_Grid_B[6][4][BN]
                             - dd2Density_Grid_B[6][4][BN] - ddDensity_Grid_B[6][2][BN];
       val[0][1] = val[0][1] - dd2Density_Grid_B[6][5][BN] - ddDensity_Grid_B[6][5][BN]
                             + dd2Density_Grid_B[6][3][BN] + ddDensity_Grid_B[6][3][BN];
       val[1][0] = val[1][0] - dd2Density_Grid_B[5][2][BN] - ddDensity_Grid_B[5][4][BN]
                             + dd2Density_Grid_B[5][4][BN] + ddDensity_Grid_B[5][2][BN];
       val[1][2] = val[1][2] + dd2Density_Grid_B[2][5][BN] + ddDensity_Grid_B[2][5][BN]
                             - dd2Density_Grid_B[2][3][BN] - ddDensity_Grid_B[2][3][BN]
                             + dd2Density_Grid_B[4][4][BN] + ddDensity_Grid_B[4][4][BN]
                             - dd2Density_Grid_B[4][2][BN] - ddDensity_Grid_B[4][2][BN];
       val[2][1] = val[2][1] - dd2Density_Grid_B[3][5][BN] - ddDensity_Grid_B[3][5][BN]
                             + dd2Density_Grid_B[3][3][BN] + ddDensity_Grid_B[3][3][BN];
       val[2][0] = val[2][0] - dd2Density_Grid_B[3][4][BN] - ddDensity_Grid_B[3][4][BN]
                             + dd2Density_Grid_B[3][2][BN] + ddDensity_Grid_B[3][2][BN];
       val[0][2] = val[0][2] + dd2Density_Grid_B[4][3][BN] + ddDensity_Grid_B[4][5][BN]
                             - dd2Density_Grid_B[4][5][BN] - ddDensity_Grid_B[4][3][BN]
                             + dd2Density_Grid_B[1][4][BN] + ddDensity_Grid_B[1][4][BN]
                             - dd2Density_Grid_B[1][2][BN] - ddDensity_Grid_B[1][2][BN];
       /* for NEGF */
       val[0][0] = val[0][0] - dd2Density_Grid_B_i[4][0][BN] + ddDensity_Grid_B_i[4][0][BN]
                             + dd2Density_Grid_B_i[4][1][BN] - ddDensity_Grid_B_i[4][1][BN];
       val[1][1] = val[1][1] - dd2Density_Grid_B_i[4][0][BN] - ddDensity_Grid_B_i[4][0][BN]
                             + dd2Density_Grid_B_i[4][1][BN] + ddDensity_Grid_B_i[4][1][BN];
       val[0][1] = val[0][1] - dd2Density_Grid_B_i[1][0][BN] - ddDensity_Grid_B_i[1][0][BN]
                             + dd2Density_Grid_B_i[1][1][BN] + ddDensity_Grid_B_i[1][1][BN];
       val[1][0] = val[1][0] + dd2Density_Grid_B_i[2][0][BN] + ddDensity_Grid_B_i[2][0][BN]
                             - dd2Density_Grid_B_i[2][1][BN] - ddDensity_Grid_B_i[2][1][BN];
       val[2][1] = val[2][1] + dd2Density_Grid_B_i[6][0][BN] - ddDensity_Grid_B_i[6][0][BN]
                             - dd2Density_Grid_B_i[6][1][BN] + ddDensity_Grid_B_i[6][1][BN];
       val[2][0] = val[2][0] + dd2Density_Grid_B_i[5][0][BN] + ddDensity_Grid_B_i[5][0][BN]
                             - dd2Density_Grid_B_i[5][1][BN] - ddDensity_Grid_B_i[5][1][BN];

       /* coef -0.5 */
       for(i=0;i<3;i++){
         for(j=0;j<3;j++){
           val[i][j] = -0.5*val[i][j];
         }
       }
    }
}

void stress(double *out, int BN)
{
    double val[3][3];
    compute_stress_full(val, BN);
    store_tensor(out, val);
}

void stress_S(double *out, int BN)
{
    double val[3][3];
    int i,j;
    compute_stress_full(val,BN);
    for (i=0; i<3; ++i) {
      for (j=i+1; j<3; ++j) {
        const double symmetric=0.5*(val[i][j]+val[j][i]);
        val[i][j]=symmetric;
        val[j][i]=symmetric;
      }
    }
    store_tensor(out, val);
}

void stress_A(double *out, int BN)
{
    double val[3][3];
    int i,j;
    compute_stress_full(val,BN);
    for (i=0; i<3; ++i) {
      val[i][i]=0.0;
      for (j=i+1; j<3; ++j) {
        const double antisymmetric=0.5*(val[i][j]-val[j][i]);
        val[i][j]=antisymmetric;
        val[j][i]=-antisymmetric;
      }
    }
    store_tensor(out, val);
}

void stress_diag(double *out, int BN)
{
    int i,j;
    double val[3][3], *val2[3], ko[3];
    compute_stress_full(val,BN);
    for (i=0; i<3; ++i) {
      for (j=i+1; j<3; ++j) {
        const double symmetric=0.5*(val[i][j]+val[j][i]);
        val[i][j]=symmetric;
        val[j][i]=symmetric;
      }
      val2[i]=val[i];
    }
    Eigen_lapack_d(val2,ko,3,3);
    store_tensor_diag(out, val, ko);
}

void dd_electron_density(double *out, int BN)
{ 
    double val[3][3];

    val[0][0] =  dd2Density_Grid_B[1][0][BN] + ddDensity_Grid_B[1][0][BN]
               + dd2Density_Grid_B[1][1][BN] + ddDensity_Grid_B[1][1][BN];
    val[1][1] =  dd2Density_Grid_B[2][0][BN] + ddDensity_Grid_B[2][0][BN]
               + dd2Density_Grid_B[2][1][BN] + ddDensity_Grid_B[2][1][BN];
    val[2][2] =  dd2Density_Grid_B[3][0][BN] + ddDensity_Grid_B[3][0][BN]
               + dd2Density_Grid_B[3][1][BN] + ddDensity_Grid_B[3][1][BN];
    val[0][1] =  dd2Density_Grid_B[4][0][BN] + ddDensity_Grid_B[4][0][BN]
               + dd2Density_Grid_B[4][1][BN] + ddDensity_Grid_B[4][1][BN];
    val[1][2] =  dd2Density_Grid_B[5][0][BN] + ddDensity_Grid_B[5][0][BN]
               + dd2Density_Grid_B[5][1][BN] + ddDensity_Grid_B[5][1][BN];
    val[2][0] =  dd2Density_Grid_B[6][0][BN] + ddDensity_Grid_B[6][0][BN]
               + dd2Density_Grid_B[6][1][BN] + ddDensity_Grid_B[6][1][BN];
    val[1][0] = val[0][1];
    val[2][1] = val[1][2];
    val[0][2] = val[2][0];

    store_tensor(out, val);
}




void stress_nonrel(double *out, int BN)   /* added by yoshida 7.January.2020 */
{ 
    int i,j;
    double val[3][3];

    if (SpinP_switch==0) {
       val[0][0] = -( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN] );
       val[1][1] = -( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN] );
       val[2][2] = -( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN] );
       val[0][1] = -( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==1 || SpinP_switch==3) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN]
                        + dd2Density_Grid_B[1][1][BN] - ddDensity_Grid_B[1][1][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN]
                        + dd2Density_Grid_B[2][1][BN] - ddDensity_Grid_B[2][1][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN]
                        + dd2Density_Grid_B[3][1][BN] - ddDensity_Grid_B[3][1][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                        + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                        + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                        + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
 /*    
    else if (SpinP_switch==3) {
       val[0][0] =   dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN]
                   + dd2Density_Grid_B[1][1][BN] - ddDensity_Grid_B[1][1][BN] ;
       val[1][1] =   dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN]
                   + dd2Density_Grid_B[2][1][BN] - ddDensity_Grid_B[2][1][BN] ;
       val[2][2] =   dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN]
                   + dd2Density_Grid_B[3][1][BN] - ddDensity_Grid_B[3][1][BN] ;
       val[0][1] =   dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN]
                   + dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] ;
       val[1][2] =   dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN]
                   + dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] ;
       val[2][0] =   dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN]
                   + dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] ;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];


       // symmetrize 
       val[0][1] = (val[0][1] + val[1][0])*0.5;
       val[1][2] = (val[1][2] + val[2][1])*0.5;
       val[2][0] = (val[2][0] + val[0][2])*0.5;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

       // coef -0.5 
       for(i=0;i<3;i++){
           for(j=0;j<3;j++){
               val[i][j] = -0.5*val[i][j];
           }
       }
 
    }
 */
    store_tensor(out, val);
}

void stress_nonrel_alpha(double *out, int BN)   /* added by yoshida 7.January.2020 */
{ 
    int i,j;
    double val[3][3];

    if (SpinP_switch==0) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==1 || SpinP_switch==3) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    store_tensor(out, val);
}

void stress_nonrel_beta(double *out, int BN)   /* added by yoshida 7.January.2020 */
{ 
    int i,j;
    double val[3][3];

    if (SpinP_switch==0) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][0][BN] - ddDensity_Grid_B[1][0][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][0][BN] - ddDensity_Grid_B[2][0][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][0][BN] - ddDensity_Grid_B[3][0][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][0][BN] - ddDensity_Grid_B[4][0][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][0][BN] - ddDensity_Grid_B[5][0][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][0][BN] - ddDensity_Grid_B[6][0][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    else if (SpinP_switch==1 || SpinP_switch==3) {
       val[0][0] = -0.5*( dd2Density_Grid_B[1][1][BN] - ddDensity_Grid_B[1][1][BN] );
       val[1][1] = -0.5*( dd2Density_Grid_B[2][1][BN] - ddDensity_Grid_B[2][1][BN] );
       val[2][2] = -0.5*( dd2Density_Grid_B[3][1][BN] - ddDensity_Grid_B[3][1][BN] );
       val[0][1] = -0.5*( dd2Density_Grid_B[4][1][BN] - ddDensity_Grid_B[4][1][BN] );
       val[1][2] = -0.5*( dd2Density_Grid_B[5][1][BN] - ddDensity_Grid_B[5][1][BN] );
       val[2][0] = -0.5*( dd2Density_Grid_B[6][1][BN] - ddDensity_Grid_B[6][1][BN] );
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
    }
    store_tensor(out, val);
}



void Eigen_lapack_d(double **a, double *ko, int n0, int EVmax)
{

  /* 
    F77_NAME(dsyevd,DSYEVD)()
  
    input:  n;
    input:  a[n][n];  matrix A
    output: a[n][n];  eigevectors
    output: ko[n];    eigenvalues 
  */
    
  static char *name="Eigen_lapack_d";

  char  *JOBZ="V";
  char  *UPLO="L";

  INTEGER n=n0;
  INTEGER LDA=n;
  double VL,VU; /* dummy */
  INTEGER IL,IU; 
  double ABSTOL=LAPACK_ABSTOL;
  INTEGER M;

  double *A;
  INTEGER LDZ=n;
  INTEGER LWORK,LIWORK;
  double *WORK;
  INTEGER *IWORK;
  INTEGER INFO;

  int i,j;

  A=(double*)malloc(sizeof(double)*n*n);

  LWORK=  1 + 6*n + 2*n*n;
  WORK=(double*)malloc(sizeof(double)*LWORK);

  LIWORK = 3 + 5*n;
  IWORK=(INTEGER*)malloc(sizeof(INTEGER)*LIWORK);


  IL = 1;
  IU = EVmax; 
 
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       //A[i*n+j]= a[i+1][j+1]; /* fukuda */
       A[i*n+j]= a[i][j];
    }
  }

#if 0
  printf("A=\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       printf("%f ",A[i*n+j]);
    }
    printf("\n");
  }
  fflush(stdout);
#endif

  F77_NAME(dsyevd,DSYEVD)( JOBZ, UPLO, &n, A, &LDA, ko, WORK, &LWORK, IWORK, &LIWORK, &INFO ); 

  /* store eigenvectors */
  for (i=0;i<EVmax;i++) {
    for (j=0;j<n;j++) {
      /* a[i+1][j+1]= Z[i*n+j]; */
      //a[j+1][i+1]= A[i*n+j]; /* fukuda */
      a[j][i]= A[i*n+j];
    }
  }

  /* comment outed by fukuda */
  ///* shift ko by 1 */
  //for (i=EVmax; i>=1; i--){
  //  ko[i]= ko[i-1];
  //}

  /*
  if (INFO>0) {
     printf("\n%s: error in dsyevd_, info=%d\n\n",name,INFO);
  }
  */

  if (INFO<0) {
     printf("%s: info=%d\n",name,INFO);
     exit(10);
  }
   
  free(IWORK); free(WORK); free(A);
}

