
/**********************************************************************
  Free_Arrays_in_Set_Grid.c:

     Free_Arrays_in_Set_Grid.c is a subrutine to free arrays.

  Log of Free_Arrays_in_Set_Grid.c:

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

void Free_Arrays_in_Set_Grid()
{

  int i,j,k;
  int spin,Mc_AN,Gc_AN,Cwan,tno0;
  int h_AN,Gh_AN,Mh_AN,Hwan,tno1,Rnh;
  int Nc,num;

  int numprocs,myid,tag=999,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************/
  /* allocated in Set_Grid.c */
  /****************************************/

    if (SpinP_switch==3){ /* spin non-collinear */
      for (k=0; k<=3; k++){
        free(Density_Grid_B[k]); 
      }
      free(Density_Grid_B); 

      for (k=0; k<3; k++){
         for (j=0; j<=5; j++){
           free(dDensity_Grid_B[k][j]); 
         }
         free(dDensity_Grid_B[k]); 
      }
      free(dDensity_Grid_B); 

      for (k=0; k<7; k++){
         for (j=0; j<=5; j++){
           free(ddDensity_Grid_B[k][j]); 
         }
         free(ddDensity_Grid_B[k]); 
      }
      free(ddDensity_Grid_B); 

      for (k=0; k<7; k++){
         for (j=0; j<=5; j++){
           free(dd2Density_Grid_B[k][j]); 
         }
         free(dd2Density_Grid_B[k]); 
      }
      free(dd2Density_Grid_B); 
    }
    else{
      for (k=0; k<=1; k++){
        free(Density_Grid_B[k]); 
      }
      free(Density_Grid_B); 

      for (k=0; k<3; k++){
         for (j=0; j<=1; j++){
           free(dDensity_Grid_B[k][j]); 
         }
         free(dDensity_Grid_B[k]); 
      }
      free(dDensity_Grid_B); 

      for (k=0; k<7; k++){
         for (j=0; j<=1; j++){
           free(ddDensity_Grid_B[k][j]); 
         }
         free(ddDensity_Grid_B[k]); 
      }
      free(ddDensity_Grid_B); 

      for (k=0; k<7; k++){
         for (j=0; j<=1; j++){
           free(dd2Density_Grid_B[k][j]); 
         }
         free(dd2Density_Grid_B[k]); 
      }
      free(dd2Density_Grid_B); 
    }

      /*** for NEGF ***/
      for (k=0; k<=1; k++){
        free(Density_Grid_B_i[k]); 
      }
      free(Density_Grid_B_i); 

      for (k=0; k<3; k++){
         for (j=0; j<=1; j++){
           free(dDensity_Grid_B_i[k][j]); 
         }
         free(dDensity_Grid_B_i[k]); 
      }
      free(dDensity_Grid_B_i); 

      for (k=0; k<7; k++){
         for (j=0; j<=1; j++){
           free(ddDensity_Grid_B_i[k][j]); 
         }
         free(ddDensity_Grid_B_i[k]); 
      }
      free(ddDensity_Grid_B_i); 

      for (k=0; k<7; k++){
         for (j=0; j<=1; j++){
           free(dd2Density_Grid_B_i[k][j]); 
         }
         free(dd2Density_Grid_B_i[k]); 
      }
      free(dd2Density_Grid_B_i); 
      /* end for NEGF */

    /* Orbs_Grid */
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      /* AITUNE */
      int Nc;
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        free(Orbs_Grid[Mc_AN][Nc]); 
      }
      free(Orbs_Grid[Mc_AN]); 
      /* AITUNE */
    }
    free(Orbs_Grid[0][0]); 
    free(Orbs_Grid[0]); 
    free(Orbs_Grid); 

    for (k=0; k<3; k++){
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        /* AITUNE */
        int Nc;
        for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
          free(dOrbs_Grid[k][Mc_AN][Nc]); 
        }
        free(dOrbs_Grid[k][Mc_AN]); 
        /* AITUNE */
      }
      free(dOrbs_Grid[k][0][0]); 
      free(dOrbs_Grid[k][0]); 
      free(dOrbs_Grid[k]); 
    }
    free(dOrbs_Grid); 

    for (k=0; k<7; k++){
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        /* AITUNE */
        int Nc;
        for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
          free(ddOrbs_Grid[k][Mc_AN][Nc]); 
        }
        free(ddOrbs_Grid[k][Mc_AN]); 
        /* AITUNE */
      }
      free(ddOrbs_Grid[k][0][0]); 
      free(ddOrbs_Grid[k][0]); 
      free(ddOrbs_Grid[k]); 
    }
    free(ddOrbs_Grid); 

    /* COrbs_Grid */
    if (Cnt_switch!=0){
      for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        for (i=0; i<Spe_Total_CNO[Cwan]; i++){
          free(COrbs_Grid[Mc_AN][i]); 
        }
        free(COrbs_Grid[Mc_AN]); 
      }
      free(COrbs_Grid[0][0]); 
      free(COrbs_Grid[0]); 
      free(COrbs_Grid); 
    }

    if (Cnt_switch!=0){
      for (k=0; k<3; k++){
        for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          for (i=0; i<Spe_Total_CNO[Cwan]; i++){
            free(CdOrbs_Grid[k][Mc_AN][i]); 
          }
          free(CdOrbs_Grid[k][Mc_AN]); 
        }
        free(CdOrbs_Grid[k][0][0]); 
        free(CdOrbs_Grid[k][0]); 
        free(CdOrbs_Grid[k]); 
      }
      free(CdOrbs_Grid); 
    }

    if (Cnt_switch!=0){
      for (k=0; k<7; k++){
        for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          for (i=0; i<Spe_Total_CNO[Cwan]; i++){
            free(CddOrbs_Grid[k][Mc_AN][i]); 
          }
          free(CddOrbs_Grid[k][Mc_AN]); 
        }
        free(CddOrbs_Grid[k][0][0]); 
        free(CddOrbs_Grid[k][0]); 
        free(CddOrbs_Grid[k]); 
      }
      free(CddOrbs_Grid); 
    }

  /* Orbs_Grid_FNAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];    

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      Gh_AN = natn[Gc_AN][h_AN];

      if (G2ID[Gh_AN]!=myid){

        Mh_AN = F_G2M[Gh_AN];
        Rnh = ncn[Gc_AN][h_AN];
        Hwan = WhatSpecies[Gh_AN];
        if (0<NumOLG[Mc_AN][h_AN]){
          int Nc;
          for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
                free(Orbs_Grid_FNAN[Mc_AN][h_AN][Nc]); 
          }
          free(Orbs_Grid_FNAN[Mc_AN][h_AN]); 
        }
      }
      else {
        free(Orbs_Grid_FNAN[Mc_AN][h_AN][0]); 
        free(Orbs_Grid_FNAN[Mc_AN][h_AN]); 
      }
    }
    free(Orbs_Grid_FNAN[Mc_AN]); 
  }
  free(Orbs_Grid_FNAN[0][0][0]); 
  free(Orbs_Grid_FNAN[0][0]); 
  free(Orbs_Grid_FNAN[0]); 
  free(Orbs_Grid_FNAN); 

  for (k=0; k<3; k++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];    

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

          Mh_AN = F_G2M[Gh_AN];
          Rnh = ncn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          if (0<NumOLG[Mc_AN][h_AN]){
            int Nc;
            for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
                  free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc]); 
            }
            free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
          }
        }
        else {
          free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN][0]); 
          free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
        }
      }
      free(dOrbs_Grid_FNAN[k][Mc_AN]); 
    }
    free(dOrbs_Grid_FNAN[k][0][0][0]); 
    free(dOrbs_Grid_FNAN[k][0][0]); 
    free(dOrbs_Grid_FNAN[k][0]); 
    free(dOrbs_Grid_FNAN[k]); 
  } /* k */
  free(dOrbs_Grid_FNAN); 

  for (k=0; k<7; k++){
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];    

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

          Mh_AN = F_G2M[Gh_AN];
          Rnh = ncn[Gc_AN][h_AN];
          Hwan = WhatSpecies[Gh_AN];
          if (0<NumOLG[Mc_AN][h_AN]){
            int Nc;
            for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
                  free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc]); 
            }
            free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
          }
        }
        else {
          free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][0]); 
          free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
        }
      }
      free(ddOrbs_Grid_FNAN[k][Mc_AN]); 
    }
    free(ddOrbs_Grid_FNAN[k][0][0][0]); 
    free(ddOrbs_Grid_FNAN[k][0][0]); 
    free(ddOrbs_Grid_FNAN[k][0]); 
    free(ddOrbs_Grid_FNAN[k]); 
  } /* k */
  free(ddOrbs_Grid_FNAN); 


  FNAN[0] = 0; 
  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

    if (Mc_AN==0) Gc_AN = 0;
    else          Gc_AN = M2G[Mc_AN];

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      free(GListTAtoms2[Mc_AN][h_AN]);
      free(GListTAtoms1[Mc_AN][h_AN]);

    }
    free(GListTAtoms2[Mc_AN]);
    free(GListTAtoms1[Mc_AN]);
  }
  free(GListTAtoms2);
  free(GListTAtoms1);


  for (ID=0; ID<numprocs; ID++){
    free(Index_Snd_Grid_A2B[ID]);
  }  
  free(Index_Snd_Grid_A2B);
  
  for (ID=0; ID<numprocs; ID++){
    free(Index_Rcv_Grid_A2B[ID]);
  }  
  free(Index_Rcv_Grid_A2B);

}
