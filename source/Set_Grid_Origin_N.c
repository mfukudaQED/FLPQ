/**********************************************************************
  Set_Grid_Origin_N.c:

     Set_Grid_Origin_N.c is a subrutine to separate lattice in the x direction.

  Log of Set_Grid_Origin_N.c:

     12/Jan/2023  Released by M.Fukuda
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"
#include "mpi.h"
#include <omp.h>


void Set_Grid_Origin_N()
{

  int i,j,k,n1;
  int ct_AN;
  double xc,yc,zc,xm,ym,zm;
  double sn1,sn2,sn3;

  if (Dens_Ngrid[0]==0 && Dens_Ngrid[1]==0 && Dens_Ngrid[2]==0){
     Dens_Ngrid[0] = Ngrid1;
     Dens_Ngrid[1] = Ngrid2;
     Dens_Ngrid[2] = Ngrid3;    
  }
  Ngrid1 = Dens_Ngrid[0];
  Ngrid2 = Dens_Ngrid[1];
  Ngrid3 = Dens_Ngrid[2];

  for (i=1; i<=3; i++){
    for (j=1; j<=3; j++){
      tv_ori[i][j] = tv[i][j];
    }
  }

  /*** lattice vectors tv are replaced by grid_vec ***/
  if(flag_grid_vec==1){
    for (i=1; i<=3; i++){
      for (j=1; j<=3; j++){
        tv[i][j] = grid_vec[i][j];
      }
    }
  }

  /* the center of the system */

  xc = 0.0;
  yc = 0.0;
  zc = 0.0;

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    xc += Gxyz[ct_AN][1];
    yc += Gxyz[ct_AN][2];
    zc += Gxyz[ct_AN][3];
  }

  xc = xc/(double)atomnum;
  yc = yc/(double)atomnum;
  zc = zc/(double)atomnum;

    /************
     start calc.
    ************/

    /* gtv */

    gtv[1][1] = tv[1][1]/(double)Ngrid1;
    gtv[1][2] = tv[1][2]/(double)Ngrid1;
    gtv[1][3] = tv[1][3]/(double)Ngrid1;

    gtv[2][1] = tv[2][1]/(double)Ngrid2;
    gtv[2][2] = tv[2][2]/(double)Ngrid2;
    gtv[2][3] = tv[2][3]/(double)Ngrid2;

    gtv[3][1] = tv[3][1]/(double)Ngrid3;
    gtv[3][2] = tv[3][2]/(double)Ngrid3;
    gtv[3][3] = tv[3][3]/(double)Ngrid3;

    sn1 = 0.5*( (Ngrid1+1) % 2 ); 
    sn2 = 0.5*( (Ngrid2+1) % 2 ); 
    sn3 = 0.5*( (Ngrid3+1) % 2 ); 

    xm = ( (double)(Ngrid1/2) - sn1 )*gtv[1][1]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][1]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][1];

    ym = ( (double)(Ngrid1/2) - sn1 )*gtv[1][2]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][2]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][2];

    zm = ( (double)(Ngrid1/2) - sn1 )*gtv[1][3]
       + ( (double)(Ngrid2/2) - sn2 )*gtv[2][3]
       + ( (double)(Ngrid3/2) - sn3 )*gtv[3][3];

    
    if ( 1.0e+8<scf_fixed_dens_origin[0] ){
      Grid_Origin[1] = xc - xm;
    }
    else{
         Grid_Origin[1] = scf_fixed_dens_origin[0];
    }
    if ( 1.0e+8<scf_fixed_dens_origin[1] ){
      Grid_Origin[2] = yc - ym;
    }
    else{
         Grid_Origin[2] = scf_fixed_dens_origin[1];
    }
    if ( 1.0e+8<scf_fixed_dens_origin[2] ){
      Grid_Origin[3] = zc - zm;
    }
    else{
         Grid_Origin[3] = scf_fixed_dens_origin[2];
    }


    Grid_Origin_N = (double**)malloc(sizeof(double*)*(Ngrid1));
    for (n1=0; n1<Ngrid1; n1++){
      Grid_Origin_N[n1] = (double*)malloc(sizeof(double)*(4));
    }

    for (n1=0; n1<Ngrid1; n1++){
      Grid_Origin_N[n1][1] = (double)n1*gtv[1][1] + Grid_Origin[1];
      Grid_Origin_N[n1][2] = (double)n1*gtv[1][2] + Grid_Origin[2];
      Grid_Origin_N[n1][3] = (double)n1*gtv[1][3] + Grid_Origin[3];
    }

}

