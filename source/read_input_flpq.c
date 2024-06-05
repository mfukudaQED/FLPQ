/**********************************************************************
  read_input_flpq.c:

     read_input_flpq.c is a subroutine to read an input file from
     flpq.inp. 

  Log of read_input_flpq.c:

     28/Jun/2016  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"

void read_input_flpq()
{
  char fname[500];
  FILE *fp;
  int i;


  flag_DMmu = 0;
  flag_DMkmu = 0;
  flag_grid_vec = 0;
    

  /****************************/
  /* input data from flpq.inp */
  /****************************/
    sprintf(fname,"flpq.inp");

    if ((fp = fopen(fname,"r")) == NULL){
      printf("cannot open %s \n", fname);
      exit(0);
    }

    
    fscanf(fp,"%lf %lf %lf \n",&scf_fixed_dens_origin[0],&scf_fixed_dens_origin[1],&scf_fixed_dens_origin[2]);
    fscanf(fp,"%d %d %d \n",&Dens_Ngrid1,&Dens_Ngrid2,&Dens_Ngrid3);

    if(flag_grid_vec==1){
      for (i=1; i<=3; i++){
        fscanf(fp,"%lf %lf %lf \n",&grid_vec[i][1],&grid_vec[i][2],&grid_vec[i][3]);
      }
    }

   fclose(fp);

}


