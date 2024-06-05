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
static void electron_density(FILE *fp, int BN, int flag_cube);
static void charge_density(FILE *fp, int BN, int flag_cube);
static void zeta_potential(FILE *fp, int BN, int flag_cube);
static void kinetic_energy_density(FILE *fp, int BN, int flag_cube);
//static void energy_density(FILE *fp, int BN, int flag_cube);
static void charge_density_alpha(FILE *fp, int BN, int flag_cube);   /* added by yoshida 23.December.2019 */
static void charge_density_beta(FILE *fp, int BN, int flag_cube);   /* added by yoshida 23.December.2019 */
static void energy(FILE *fp, int BN, int flag_cube);   /* added by takamatsu 19.February.2020 */
static void chemical_potential(FILE *fp, int BN, int flag_cube);   /* added by takamatsu 19.February.2020 */

/*  vector  */
static void charge_current(FILE *fp, int BN);
static void kinetic_momentum(FILE *fp, int BN);
static void kinetic_momentum_alpha(FILE *fp, int BN);
static void kinetic_momentum_beta(FILE *fp, int BN);
//static void kinetic_momentum_spinup(FILE *fp, int BN);
//static void kinetic_momentum_spindown(FILE *fp, int BN);
static void spin_vorticity(FILE *fp, int BN);
static void spin_angular_momentum(FILE *fp, int BN);
//static void spin_torque(FILE *fp, int BN);
static void zeta_force(FILE *fp, int BN);
static void pau_dd(FILE *fp, int BN);
//static void tension(FILE *fp, int BN);
//static void lorentz_force(FILE *fp, int BN);
static void gradient_electron_density(FILE *fp, int BN);

static void spin_z_current(FILE *fp, int BN);   /* added by kuroda 7.June.2019 */
static void charge_current_alpha(FILE *fp, int BN);     /* added by kuroda 5.December.2019 */
static void charge_current_beta(FILE *fp, int BN);   /* added by kuroda 5.December.2019 */
static void velocity(FILE *fp, int BN);   /* added by kuroda 31.January.2020 */
static void torq_zeta(FILE *fp, int BN);      /* added by shimada 27.December.2019 */
static void soc_part1(FILE *fp, int BN);       /* added by shimada September.2020 */
static void soc_part2(FILE *fp, int BN);       /* added by shimada September.2020 */

/*  tensor  */
static void stress(FILE *fp, int BN);
static void stress_S(FILE *fp, int BN);
static void stress_A(FILE *fp, int BN);
static void stress_diag(FILE *fp, int BN);
static void dd_electron_density(FILE *fp, int BN);
static void stress_nonrel(FILE *fp, int BN);   /* added by yoshida 7.January.2020 */
static void stress_nonrel_alpha(FILE *fp, int BN);   /* added by yoshida 7.January.2020 */
static void stress_nonrel_beta(FILE *fp, int BN);   /* added by yoshida 7.January.2020 */

/* diagonalization */
static void Eigen_lapack_d(double **a, double *ko, int n, int EVmax);

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

void scalar_density(int n, FILE *fp, int BN, int flag_cube)
{
  if (n==0) {
     electron_density(fp,BN,flag_cube);
  }
  else if (n==1) {
     charge_density(fp,BN,flag_cube);
  }
  else if (n==2) {
     zeta_potential(fp,BN,flag_cube);
  }
  else if (n==3) {
     kinetic_energy_density(fp,BN,flag_cube);
  }
  //else if (n==4) {
  //   energy_density(fp,BN,flag_cube);
  //}
  else if (n==5) {
     charge_density_alpha(fp,BN,flag_cube);  /* added by yoshida 23.December.2019 */
  }
  else if (n==6) {
     charge_density_beta(fp,BN,flag_cube);  /* added by yoshida 23.December.2019 */
  }
  else if (n==7) {
     energy(fp,BN,flag_cube);  /* added by takamatsu 19.February.2020 */
  }
  else if (n==8) {
     chemical_potential(fp,BN,flag_cube);  /* added by takamatsu 19.February.2020 */
  }
}

void vector_density(int n, FILE *fp, int BN)
{
  if (n==0) {
     charge_current(fp,BN);
  }
  else if (n==1) {
     kinetic_momentum(fp,BN);
  }
  else if (n==2) {
     spin_vorticity(fp,BN);
  }
  else if (n==3) {
     spin_angular_momentum(fp,BN);
  }
  //else if (n==4) {
  //   spin_torque(fp,BN);
  //}
  else if (n==5) {
     zeta_force(fp,BN);
  }
  else if (n==6) {
     pau_dd(fp,BN);
  }
  //else if (n==7) {
  //   tension(fp,BN);
  //}
  //else if (n==8) {
  //   lorentz_force(fp,BN);
  //}
  else if (n==9) {
     gradient_electron_density(fp,BN);
  }
  
    else if (n==10) {
     spin_z_current(fp,BN);  /* added by kuroda 7.June.2019 */
  }
  
    else if (n==11) {
     charge_current_alpha(fp,BN);  /* added by kuroda 5.December.2019 */
  }
  
    else if (n==12) {
     charge_current_beta(fp,BN);  /* added by kuroda 5.December.2019 */
  }

    else if (n==13) {
     velocity(fp,BN);  /* added by kuroda 31.January.2020 */
  }

    else if (n==14) {
     torq_zeta(fp,BN);  /* added by shimada 27.December.2019 */
  }

    else if (n==15) {
     soc_part1(fp,BN);  /* added by shimada September2020 */
  }

    else if (n==16) {
     kinetic_momentum_alpha(fp,BN);
  }
    else if (n==17) {
     kinetic_momentum_beta(fp,BN);
  }
   else if (n==18) {
     soc_part2(fp,BN);  /* added by shimada September2020 */
  }
}

void tensor_density(int n, FILE *fp, int BN)
{
  if (n==0) {
     stress(fp,BN);
  }
  else if (n==1) {
     stress_S(fp,BN);
  }
  else if (n==2) {
     stress_A(fp,BN);
  }
  else if (n==3) {
     stress_diag(fp,BN);
  }
  else if (n==4) {
     dd_electron_density(fp,BN);
  }
  else if (n==5) {
     stress_nonrel(fp,BN);   /* added by yoshida 7.January.2020 */
  }
  else if (n==6) {
     stress_nonrel_alpha(fp,BN);   /* added by yoshida 7.January.2020 */
  }
  else if (n==7) {
     stress_nonrel_beta(fp,BN);   /* added by yoshida 7.January.2020 */
  }
}

/******************************************
        Local physical quantities
 ******************************************/
/*******************************
             scalar
 *******************************/
void electron_density(FILE *fp, int BN, int flag_cube)
{ 
    double val;

    if (SpinP_switch==0) {
       val = 2.0*Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN] + Density_Grid_B[1][BN];
    }
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);

}

void charge_density(FILE *fp, int BN, int flag_cube)
{ 
    double val;

    if (SpinP_switch==0) {
       val = 2.0*Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN] + Density_Grid_B[1][BN];
    }
    val = - val;
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);

}

void zeta_potential(FILE *fp, int BN, int flag_cube)
{ 
    double val;

    val =  dDensity_Grid_B[0][3][BN] - dDensity_Grid_B[0][5][BN]
          -dDensity_Grid_B[1][2][BN] + dDensity_Grid_B[1][4][BN]
          +dDensity_Grid_B_i[2][0][BN] - dDensity_Grid_B_i[2][1][BN];
    val = 0.5*val;
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);

}

void kinetic_energy_density(FILE *fp, int BN, int flag_cube)
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
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);

}


void charge_density_alpha(FILE *fp, int BN, int flag_cube)   /* added by yoshida 23.December.2019 */
{ 
    double val;

    if (SpinP_switch==0) {
       val = Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[0][BN];
    }
    val = - val;
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);
}

void charge_density_beta(FILE *fp, int BN, int flag_cube)   /* added by yoshida 23.December.2019 */
{ 
    double val;

    if (SpinP_switch==0) {
       val = Density_Grid_B[0][BN];
    }
    else {
       val = Density_Grid_B[1][BN];
    }
    val = - val;
    print_scalar(fp, val, flag_cube);
    //printf("%24.14E %d\n", val, BN);
}

void energy(FILE *fp, int BN, int flag_cube)   /* added by takamatsu 19.February.2020 */
{ 
//start of calculation of strss_diag----------------------------------------------------------------------
    int i,j;
    double val[3][3];
    double *val2[3];
    double ko[3];

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
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
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

       /* symmetrize */
       val[0][1] = (val[0][1] + val[1][0])*0.5;
       val[1][2] = (val[1][2] + val[2][1])*0.5;
       val[2][0] = (val[2][0] + val[0][2])*0.5;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

    /* coef -0.5 */
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            val[i][j] = -0.5*val[i][j];
        }
    }
    }

    /* diagonalize */
    for(i=0;i<3;i++) val2[i] = val[i]; /* set pointer */
    Eigen_lapack_d(val2,ko,3,3);
//end of calculation of strss_diag----------------------------------------------------------------------

    double energy = 0.5*(ko[0] + ko[1] + ko[2]);

    print_scalar(fp, energy, flag_cube);
    //printf("%24.14E %d\n", val, BN);
}

void chemical_potential(FILE *fp, int BN, int flag_cube)   /* added by takamatsu 19.February.2020 */
{ 
//start of calculation of strss_diag----------------------------------------------------------------------
    int i,j;
    double val[3][3];
    double *val2[3];
    double ko[3];

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
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
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

       /* symmetrize */
       val[0][1] = (val[0][1] + val[1][0])*0.5;
       val[1][2] = (val[1][2] + val[2][1])*0.5;
       val[2][0] = (val[2][0] + val[0][2])*0.5;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

    /* coef -0.5 */
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            val[i][j] = -0.5*val[i][j];
        }
    }
    }

    /* diagonalize */
    for(i=0;i<3;i++) val2[i] = val[i]; /* set pointer */
    Eigen_lapack_d(val2,ko,3,3);
//end of calculation of strss_diag----------------------------------------------------------------------

    double energy = 0.5*(ko[0] + ko[1] + ko[2]);

//start of calculation of electron_density-----------------------------------------------------------

    double dens;

    if (SpinP_switch==0) {
       dens = 2.0*Density_Grid_B[0][BN];
    }
    else {
       dens = Density_Grid_B[0][BN] + Density_Grid_B[1][BN];
    }

//end of calculation of electron_density-----------------------------------------------------------

    double chem_pot;
    double dens_threashold;
    dens_threashold = 1.e-16;

    if(fabs(dens)<dens_threashold){
      chem_pot = 0.0;
    }
    else{
      chem_pot = energy / dens;
    }

    print_scalar(fp, chem_pot, flag_cube);
    //printf("%24.14E %d\n", val, BN);
}


/*******************************
            vector 
 *******************************/
void charge_current(FILE *fp, int BN)
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
    print_vector(fp, val);
}

void kinetic_momentum(FILE *fp, int BN)
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
    print_vector(fp, val);
}

void kinetic_momentum_alpha(FILE *fp, int BN)
{
    double val[3];
    
    val[0] = dDensity_Grid_B_i[0][0][BN];
    val[1] = dDensity_Grid_B_i[1][0][BN];
    val[2] = dDensity_Grid_B_i[2][0][BN];
    print_vector(fp, val);
}

void kinetic_momentum_beta(FILE *fp, int BN)
{
    double val[3];
    
    val[0] = dDensity_Grid_B_i[0][1][BN];
    val[1] = dDensity_Grid_B_i[1][1][BN];
    val[2] = dDensity_Grid_B_i[2][1][BN];
    print_vector(fp, val);
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

void spin_vorticity(FILE *fp, int BN)
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
    print_vector(fp, val);
}

void spin_angular_momentum(FILE *fp, int BN)
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
    print_vector(fp, val);
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

void torq_zeta(FILE *fp, int BN) /* added by matsumoto shimada 27.December.2019 */
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
    print_vector(fp, a1);
}
void zeta_force(FILE *fp, int BN)
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
    print_vector(fp, val);
}

void pau_dd(FILE *fp, int BN)
{
    double val[3];

    val[0] = -0.5*( ddDensity_Grid_B[0][3][BN] - ddDensity_Grid_B[0][5][BN] );
    val[1] = -0.5*(-ddDensity_Grid_B[0][2][BN] + ddDensity_Grid_B[0][4][BN] );
    val[2] = -0.5*(-ddDensity_Grid_B_i[0][0][BN] - ddDensity_Grid_B_i[0][1][BN] );
    print_vector(fp, val);
}

void gradient_electron_density(FILE *fp, int BN)
{
    double val[3];

    val[0] = 2.0*( dDensity_Grid_B[0][0][BN] + dDensity_Grid_B[0][1][BN] );
    val[1] = 2.0*( dDensity_Grid_B[1][0][BN] + dDensity_Grid_B[1][1][BN] );
    val[2] = 2.0*( dDensity_Grid_B[2][0][BN] + dDensity_Grid_B[2][1][BN] );
    print_vector(fp, val);
}


void spin_z_current(FILE *fp, int BN)  /* added by kuroda 7.June.2019 */
{
    double val[3];
    
    val[0] = 0.5*(dDensity_Grid_B_i[0][0][BN] - dDensity_Grid_B_i[0][1][BN]);
    val[1] = 0.5*(dDensity_Grid_B_i[1][0][BN] - dDensity_Grid_B_i[1][1][BN]);
    val[2] = 0.5*(dDensity_Grid_B_i[2][0][BN] - dDensity_Grid_B_i[2][1][BN]);
    print_vector(fp, val);
}


void charge_current_alpha(FILE *fp, int BN)  /* added by kuroda 5.December.2019 */
{
    double val[3];
    
    val[0] = -dDensity_Grid_B_i[0][0][BN];
    val[1] = -dDensity_Grid_B_i[1][0][BN];
    val[2] = -dDensity_Grid_B_i[2][0][BN];
    print_vector(fp, val);
}

void charge_current_beta(FILE *fp, int BN)  /* added by kuroda 5.December.2019 */
{
    double val[3];
    
    val[0] = -dDensity_Grid_B_i[0][1][BN];
    val[1] = -dDensity_Grid_B_i[1][1][BN];
    val[2] = -dDensity_Grid_B_i[2][1][BN];
    print_vector(fp, val);
}


void velocity(FILE *fp, int BN)  /* added by kuroda 31.January.2020 */
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

    print_vector(fp, val);
}

void soc_part1(FILE *fp, int BN) /* added by shimada September.2020 */
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
    print_vector(fp, val);
}

void soc_part2(FILE *fp, int BN) /* added by shimada September.2020 */
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
    print_vector(fp, val);
}
/*******************************
             tensor
 *******************************/
void stress(FILE *fp, int BN)
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

       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
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
    print_tensor(fp, val);
}

void stress_S(FILE *fp, int BN)
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

       /* for NEGF */
       if (Solver==4){
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

       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
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

       /* symmetrize */
       val[0][1] = (val[0][1] + val[1][0])*0.5;
       val[1][2] = (val[1][2] + val[2][1])*0.5;
       val[2][0] = (val[2][0] + val[0][2])*0.5;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

       /* coef -0.5 */
       for(i=0;i<3;i++){
         for(j=0;j<3;j++){
           val[i][j] = -0.5*val[i][j];
         }
       }
    }
    print_tensor(fp, val);
}

void stress_A(FILE *fp, int BN)
{ 
    int i,j;
    double val[3][3];

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
       val[1][0] = -val[0][1];
       val[2][1] = -val[1][2];
       val[0][2] = -val[2][0];

       /* coef -0.5 */
       for(i=0;i<3;i++){
         for(j=0;j<3;j++){
           val[i][j] = -0.5*val[i][j];
         }
       }
    }
    print_tensor(fp, val);
}

void stress_diag(FILE *fp, int BN)
{ 
    int i,j;
    double val[3][3];
    double *val2[3];
    double ko[3];

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

       if (Solver==4) {
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
       }

       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];
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

       /* symmetrize */
       val[0][1] = (val[0][1] + val[1][0])*0.5;
       val[1][2] = (val[1][2] + val[2][1])*0.5;
       val[2][0] = (val[2][0] + val[0][2])*0.5;
       val[1][0] = val[0][1];
       val[2][1] = val[1][2];
       val[0][2] = val[2][0];

       /* coef -0.5 */
       for(i=0;i<3;i++){
         for(j=0;j<3;j++){
           val[i][j] = -0.5*val[i][j];
         }
       }
    }
    /* diagonalize */
    for(i=0;i<3;i++) val2[i] = val[i]; /* set pointer */
    Eigen_lapack_d(val2,ko,3,3);
    print_tensor_diag(fp,val,ko);
}

void dd_electron_density(FILE *fp, int BN)
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

    print_tensor(fp, val);
}




void stress_nonrel(FILE *fp, int BN)   /* added by yoshida 7.January.2020 */
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
    print_tensor(fp, val);
}

void stress_nonrel_alpha(FILE *fp, int BN)   /* added by yoshida 7.January.2020 */
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
    print_tensor(fp, val);
}

void stress_nonrel_beta(FILE *fp, int BN)   /* added by yoshida 7.January.2020 */
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
    print_tensor(fp, val);
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

