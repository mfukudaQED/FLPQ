/**********************************************************************
  Free_Arrays.c:

     Free_Arrays.c is a subrutine to free arrays.

  Log of Free_Arrays.c:

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

void Free_Arrays()
{

  int i,j,k,n1;
  int spin,Mc_AN,Gc_AN,Cwan,tno0;
  int h_AN,Gh_AN,Mh_AN,Hwan,tno1,Rnh;

  int numprocs,myid,tag=999,ID;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /********************************/
  /* allocated in read_input_DM.c */
  /********************************/


  for (spin=0; spin<=SpinP_switch; spin++){
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
        Gc_AN = 0;
        tno0 = 1;
      }
      else{
        Gc_AN = M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        tno0 = Spe_Total_NO[Cwan];  
      }    

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        for (i=0; i<tno0; i++){
          free(DM[0][spin][Mc_AN][h_AN][i]); 
        }
        free(DM[0][spin][Mc_AN][h_AN]); 
      }
      free(DM[0][spin][Mc_AN]);
    }
    free(DM[0][spin]); 
  }
  free(DM[0]); 
  free(DM); 

  if ((Solver==4)||(SpinP_switch>1)){
    for (spin=0; spin<=SpinP_switch; spin++){
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

        if (Mc_AN==0){
          Gc_AN = 0;
          tno0 = 1;
        }
        else{
          Gc_AN = M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          tno0 = Spe_Total_NO[Cwan];  
        }    

        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          for (i=0; i<tno0; i++){
            free(iDM[0][spin][Mc_AN][h_AN][i]); 
          }
          free(iDM[0][spin][Mc_AN][h_AN]); 
        }
        free(iDM[0][spin][Mc_AN]);
      }
      free(iDM[0][spin]); 
    }
    free(iDM[0]); 
    free(iDM); 
  }

  ///****************************************/
  ///* allocated in Set_Grid.c */
  ///****************************************/


  //if (SpinP_switch==3){ /* spin non-collinear */
  //  for (k=0; k<=3; k++){
  //    free(Density_Grid_B[k]); 
  //  }
  //  free(Density_Grid_B); 
  //  for (k=0; k<3; k++){
  //     for (j=0; j<=5; j++){
  //       free(dDensity_Grid_B[k][j]); 
  //     }
  //     free(dDensity_Grid_B[k]); 
  //  }
  //  free(dDensity_Grid_B); 
  //  free(ddDensity_Grid_B); 
  //  for (k=0; k<7; k++){
  //     for (j=0; j<=5; j++){
  //       free(ddDensity_Grid_B[k][j]); 
  //     }
  //     free(ddDensity_Grid_B[k]); 
  //  }
  //  for (k=0; k<7; k++){
  //     for (j=0; j<=5; j++){
  //       free(dd2Density_Grid_B[k][j]); 
  //     }
  //     free(dd2Density_Grid_B[k]); 
  //  }
  //  free(dd2Density_Grid_B); 
  //}
  //else{
  //  for (k=0; k<=1; k++){
  //    free(Density_Grid_B[k]); 
  //  }
  //  free(Density_Grid_B); 
  //  for (k=0; k<3; k++){
  //     for (j=0; j<=1; j++){
  //       free(dDensity_Grid_B[k][j]); 
  //     }
  //     free(dDensity_Grid_B[k]); 
  //  }
  //  free(dDensity_Grid_B); 
  //  for (k=0; k<7; k++){
  //     for (j=0; j<=1; j++){
  //       free(ddDensity_Grid_B[k][j]); 
  //     }
  //     free(ddDensity_Grid_B[k]); 
  //  }
  //  free(ddDensity_Grid_B); 
  //  for (k=0; k<7; k++){
  //     for (j=0; j<=1; j++){
  //       free(dd2Density_Grid_B[k][j]); 
  //     }
  //     free(dd2Density_Grid_B[k]); 
  //  }
  //  free(dd2Density_Grid_B); 
  //}

  //  /*** for NEGF ***/
  //  for (k=0; k<=1; k++){
  //    free(Density_Grid_B_i[k]); 
  //  }
  //  free(Density_Grid_B_i); 
  //  for (k=0; k<3; k++){
  //     for (j=0; j<=1; j++){
  //       free(dDensity_Grid_B_i[k][j]); 
  //     }
  //     free(dDensity_Grid_B_i[k]); 
  //  }
  //  free(dDensity_Grid_B_i); 

  //  for (k=0; k<7; k++){
  //     for (j=0; j<=1; j++){
  //       free(ddDensity_Grid_B_i[k][j]); 
  //     }
  //     free(ddDensity_Grid_B_i[k]); 
  //  }
  //  free(ddDensity_Grid_B_i); 
  //  for (k=0; k<7; k++){
  //     free(dd2Density_Grid_B_i[k]); 
  //     for (j=0; j<=1; j++){
  //       free(dd2Density_Grid_B_i[k][j]); 
  //     }
  //  }
  //  free(dd2Density_Grid_B_i); 

  ///* Orbs_Grid */
  //for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //  Gc_AN = F_M2G[Mc_AN];
  //  Cwan = WhatSpecies[Gc_AN];
  //  int Nc;
  //  for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
  //    free(Orbs_Grid[Mc_AN][Nc]); 
  //  }
  //  free(Orbs_Grid[Mc_AN]); 
  //}
  //free(Orbs_Grid[0][0]); 
  //free(Orbs_Grid[0]); 
  //free(Orbs_Grid); 

  //for (k=0; k<3; k++){
  //  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //    Gc_AN = F_M2G[Mc_AN];
  //    Cwan = WhatSpecies[Gc_AN];
  //    int Nc;
  //    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
  //      free(dOrbs_Grid[k][Mc_AN][Nc]); 
  //    }
  //    free(dOrbs_Grid[k][Mc_AN]); 
  //  }
  //  free(dOrbs_Grid[k][0][0]); 
  //  free(dOrbs_Grid[k][0]); 
  //  free(dOrbs_Grid[k]); 
  //}
  //free(dOrbs_Grid); 

  //for (k=0; k<7; k++){
  //  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //    Gc_AN = F_M2G[Mc_AN];
  //    Cwan = WhatSpecies[Gc_AN];
  //    int Nc;
  //    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
  //      free(ddOrbs_Grid[k][Mc_AN][Nc]); 
  //    }
  //    free(ddOrbs_Grid[k][Mc_AN]); 
  //  }
  //  free(ddOrbs_Grid[k][0][0]); 
  //  free(ddOrbs_Grid[k][0]); 
  //  free(ddOrbs_Grid[k]); 
  //}
  //free(ddOrbs_Grid); 

  ///* COrbs_Grid */
  //if (Cnt_switch!=0){
  //  for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
  //    Gc_AN = F_M2G[Mc_AN];
  //    Cwan = WhatSpecies[Gc_AN];
  //    for (i=0; i<Spe_Total_CNO[Cwan]; i++){
  //      free(COrbs_Grid[Mc_AN][i]); 
  //    }
  //    free(COrbs_Grid[Mc_AN]); 
  //  }
  //  free(COrbs_Grid[0][0]); 
  //  free(COrbs_Grid[0]); 
  //  free(COrbs_Grid); 
  //}

  //if (Cnt_switch!=0){
  //  for (k=0; k<3; k++){
  //    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
  //      Gc_AN = F_M2G[Mc_AN];
  //      Cwan = WhatSpecies[Gc_AN];
  //      for (i=0; i<Spe_Total_CNO[Cwan]; i++){
  //        free(CdOrbs_Grid[k][Mc_AN][i]); 
  //      }
  //      free(CdOrbs_Grid[k][Mc_AN]); 
  //    }
  //    free(CdOrbs_Grid[k][0][0]); 
  //    free(CdOrbs_Grid[k][0]); 
  //    free(CdOrbs_Grid[k]); 
  //  }
  //  free(CdOrbs_Grid); 
  //}

  //if (Cnt_switch!=0){
  //  for (k=0; k<7; k++){
  //    for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
  //      Gc_AN = F_M2G[Mc_AN];
  //      Cwan = WhatSpecies[Gc_AN];
  //      free(CddOrbs_Grid[k][Mc_AN]); 
  //      for (i=0; i<Spe_Total_CNO[Cwan]; i++){
  //        free(CddOrbs_Grid[k][Mc_AN][i]); 
  //      }
  //    }
  //    free(CddOrbs_Grid[k][0][0]); 
  //    free(CddOrbs_Grid[k][0]); 
  //    free(CddOrbs_Grid[k]); 
  //  }
  //  free(CddOrbs_Grid); 
  //}

  ///* Orbs_Grid_FNAN */

  //for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //  Gc_AN = M2G[Mc_AN];    

  //  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
  //    Gh_AN = natn[Gc_AN][h_AN];

  //    if (G2ID[Gh_AN]!=myid){

  //      Mh_AN = F_G2M[Gh_AN];
  //      Rnh = ncn[Gc_AN][h_AN];
  //      Hwan = WhatSpecies[Gh_AN];
  //      if (0<NumOLG[Mc_AN][h_AN]){
  //        int Nc;
  //        for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
  //              free(Orbs_Grid_FNAN[Mc_AN][h_AN][Nc]); 
  //        }
  //        free(Orbs_Grid_FNAN[Mc_AN][h_AN]); 
  //      }
  //    }
  //    else {
  //      free(Orbs_Grid_FNAN[Mc_AN][h_AN][0]); 
  //      free(Orbs_Grid_FNAN[Mc_AN][h_AN]); 
  //    }
  //  }
  //  free(Orbs_Grid_FNAN[Mc_AN]); 
  //}
  //free(Orbs_Grid_FNAN[0][0][0]); 
  //free(Orbs_Grid_FNAN[0][0]); 
  //free(Orbs_Grid_FNAN[0]); 
  //free(Orbs_Grid_FNAN); 

  //for (k=0; k<3; k++){
  //  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //    Gc_AN = M2G[Mc_AN];    

  //    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
  //      Gh_AN = natn[Gc_AN][h_AN];

  //      if (G2ID[Gh_AN]!=myid){

  //        Mh_AN = F_G2M[Gh_AN];
  //        Rnh = ncn[Gc_AN][h_AN];
  //        Hwan = WhatSpecies[Gh_AN];
  //        if (0<NumOLG[Mc_AN][h_AN]){
  //          int Nc;
  //          for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
  //                free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc]); 
  //          }
  //          free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
  //        }
  //      }
  //      else {
  //        free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN][0]); 
  //        free(dOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
  //      }
  //    }
  //    free(dOrbs_Grid_FNAN[k][Mc_AN]); 
  //  }
  //  free(dOrbs_Grid_FNAN[k][0][0][0]); 
  //  free(dOrbs_Grid_FNAN[k][0][0]); 
  //  free(dOrbs_Grid_FNAN[k][0]); 
  //  free(dOrbs_Grid_FNAN[k]); 
  //} /* k */
  //free(dOrbs_Grid_FNAN); 

  //for (k=0; k<7; k++){
  //  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
  //    Gc_AN = M2G[Mc_AN];    

  //    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
  //      Gh_AN = natn[Gc_AN][h_AN];

  //      if (G2ID[Gh_AN]!=myid){

  //        Mh_AN = F_G2M[Gh_AN];
  //        Rnh = ncn[Gc_AN][h_AN];
  //        Hwan = WhatSpecies[Gh_AN];
  //        if (0<NumOLG[Mc_AN][h_AN]){
  //          int Nc;
  //          for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
  //                free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc]); 
  //          }
  //          free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
  //        }
  //      }
  //      else {
  //        free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][0]); 
  //        free(ddOrbs_Grid_FNAN[k][Mc_AN][h_AN]); 
  //      }
  //    }
  //    free(ddOrbs_Grid_FNAN[k][Mc_AN]); 
  //  }
  //  free(ddOrbs_Grid_FNAN[k][0][0][0]); 
  //  free(ddOrbs_Grid_FNAN[k][0][0]); 
  //  free(ddOrbs_Grid_FNAN[k][0]); 
  //  free(ddOrbs_Grid_FNAN[k]); 
  //} /* k */
  //free(ddOrbs_Grid_FNAN); 


  //FNAN[0] = 0; 
  //for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

  //  if (Mc_AN==0) Gc_AN = 0;
  //  else          Gc_AN = M2G[Mc_AN];

  //  for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
  //    free(GListTAtoms2[Mc_AN][h_AN]);
  //    free(GListTAtoms1[Mc_AN][h_AN]);

  //  }
  //  free(GListTAtoms2[Mc_AN]);
  //  free(GListTAtoms1[Mc_AN]);
  //}
  //free(GListTAtoms2);
  //free(GListTAtoms1);


  //for (ID=0; ID<numprocs; ID++){
  //  free(Index_Snd_Grid_A2B[ID]);
  //}  
  //free(Index_Snd_Grid_A2B);
  //
  //for (ID=0; ID<numprocs; ID++){
  //  free(Index_Rcv_Grid_A2B[ID]);
  //}  
  //free(Index_Rcv_Grid_A2B);




  /****************************************/
  /* allocated in read_input_openmx_lpq.c */
  /****************************************/

  for (i=0; i<(atomnum+1); i++){
    free(Gxyz[i]);
  }
  free(Gxyz);

  free(G2ID);

  for (i=0; i<SpeciesNum; i++){
    free(Spe_Num_Basis[i]);
  }  
  free(Spe_Num_Basis);
   
  for (i=0; i<SpeciesNum; i++){
    free(Spe_Num_CBasis[i]);
  }  
  free(Spe_Num_CBasis);
  
  for (i=0; i<List_YOUSO[18]; i++){
    free(Spe_PAO_RV[i]);
  }
  free(Spe_PAO_RV);

  for (i=0; i<List_YOUSO[18]; i++){
    for (j=0; j<=List_YOUSO[25]; j++){
      for (k=0; k<List_YOUSO[24]; k++){
        free(Spe_PAO_RWF[i][j][k]);
      }
      free(Spe_PAO_RWF[i][j]);
    }
    free(Spe_PAO_RWF[i]);
  }
  free(Spe_PAO_RWF);

  free(Spe_Atom_Cut1);

  for (i=0; i<SpeciesNum; i++){
    free(Spe_Specified_Num[i]); 
  }
  free(Spe_Specified_Num); 

  for (i=0; i<SpeciesNum; i++){
    for (j=0; j<Spe_Total_NO[i]; j++){
      free(Spe_Trans_Orbital[i][j]); 
    }
    free(Spe_Trans_Orbital[i]); 
  }
  free(Spe_Trans_Orbital); 

  /* for Print_CubeTitle */
  free(Spe_WhatAtom);

  /* for Print_CubeTitle */
  free(Spe_Core_Charge);

  /* for Print_CubeTitle */
  free(InitN_USpin);
  free(InitN_DSpin);
  free(WhatSpecies);
  free(GridN_Atom);
  for (i=0; i<=atomnum; i++){
    free(ncn[i]);
  }
  free(ncn);


  free(Spe_Total_NO);
  free(Spe_Total_CNO);
  free(Spe_Num_Mesh_PAO);
  free(Spe_MaxL_Basis);

  free(F_Snd_Num);
  free(F_Rcv_Num);
  free(F_TopMAN);
  free(Num_Snd_Grid_A2B);
  free(Num_Rcv_Grid_A2B);


  {
  int n,TN;
  TN = (2*CpyCell+1)*(2*CpyCell+1)*(2*CpyCell+1) - 1;

    for (i=0; i<(TN+1); i++){
      free(atv[i]);
    }
    free(atv);

    n = 2*CpyCell + 4;
    for (i=0; i<n; i++){
      for (j=0; j<n; j++){
        free(ratv[i][j]);
      }
      free(ratv[i]);
    }
    free(ratv);

    for (i=0; i<(TN+1); i++){
      free(atv_ijk[i]);
    }
    free(atv_ijk);
  }

  for (i=0; i<=atomnum; i++){
    free(natn[i]);
  }
  free(natn);


   /*  NumOLG is setted in Set_Grid.c */
  {
    int Mc_AN;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(NumOLG[Mc_AN]);
    }
    free(NumOLG);
  }


  if (Cnt_switch==1){
    for (i=0; i<=(Matomnum+MatomnumF); i++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(CntCoes[i][j]);
      }
      free(CntCoes[i]);
    }
    free(CntCoes);
  }

  free(FNAN);

  /****************************************/
  /* allocated in Set_Inf_SndRcv.c */
  /****************************************/

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_GAN[ID]);
  }
  free(Rcv_GAN);

  for (ID=0; ID<numprocs; ID++){
    free(Snd_MAN[ID]);
  }
  free(Snd_MAN);

  for (ID=0; ID<numprocs; ID++){
    free(Snd_GAN[ID]);
  }
  free(Snd_GAN);

  free(F_M2G);
  free(F_G2M);


  /****************************************/
  /* allocated in Set_Allocate_Atom2CPU.c */
  /****************************************/

  free(M2G);

  /****************************************/
  /* allocated in Set_Grid_Origin_N.c */
  /****************************************/
    for (n1=0; n1<Ngrid1; n1++){
      free(Grid_Origin_N[n1]);
    }
    free(Grid_Origin_N);

}
