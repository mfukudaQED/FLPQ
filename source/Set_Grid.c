/**********************************************************************
  Set_Grid.c:

     Set_Grid.c is a subroutine to cut grids for
     calculating local physical quantities. 

  Log of Set_Grid.c:

     28/Jun/2016  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "omp.h"
#include "flpq.h"

#define  measure_time   0

static void UCell_Box(int MD_iter, int estimate_switch, int CpyCell);
static void Construct_MPI_Data_Structure_Grid();

void Set_Grid()
{
  int My_Max_GridN_Atom,My_Max_OneD_Grids,My_Max_NumOLG;
  int numprocs,myid,tag=999,ID;
  int i,j,k,m,l,ct_AN,h_AN,Gh_AN,Mc_AN,Gc_AN,s1,s2;
  int tno,tno0,tno1,Cwan,Hwan,so,nc,ns,spin;
  int num,wan,n2,wanA,Gi,Max_Num_Cells0;
  int NO1,Mh_AN,Rnh;
  int size_Orbs_Grid,size_COrbs_Grid;
  int size_dOrbs_Grid,size_CdOrbs_Grid; /* fukuda */
  int size_ddOrbs_Grid,size_CddOrbs_Grid; /* fukuda */
  int size_Orbs_Grid_FNAN;
  int size_dOrbs_Grid_FNAN; /* fukuda */
  int size_ddOrbs_Grid_FNAN; /* fukuda */
  double stime,etime;
  double TStime,TEtime;

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  //if (myid==Host_ID && 0<level_stdout){
  //  printf("\n*******************************************************\n"); 
  //  printf("        Analysis of neigbors and setting of grids (in Reset UCell Box)  \n");
  //  printf("*******************************************************\n\n"); 
  //}

  /****************************************************
                       UCell_Box
  ****************************************************/
  //if (Dens_Ngrid[0]==0 && Dens_Ngrid[1]==0 && Dens_Ngrid[2]==0){
  //   Dens_Ngrid[0] = Ngrid1;
  //   Dens_Ngrid[1] = Ngrid2;
  //   Dens_Ngrid[2] = Ngrid3;    
  //}
  //Ngrid1 = Dens_Ngrid[0];
  //Ngrid2 = Dens_Ngrid[1];
  //Ngrid3 = Dens_Ngrid[2];

  //for (i=1; i<=3; i++){
  //  for (j=1; j<=3; j++){
  //    tv_ori[i][j] = tv[i][j];
  //  }
  //}

  ///*** lattice vectors tv are replaced by grid_vec ***/
  //if(flag_grid_vec==1){
  //  for (i=1; i<=3; i++){
  //    for (j=1; j<=3; j++){
  //      tv[i][j] = grid_vec[i][j];
  //    }
  //  }
  //}

//  if (MD_iter==1 && UCell_flag==1){

    /*************************************
      find 
            Max_GridN_Atom
            Max_OneD_Grids
    *************************************/

    Max_GridN_Atom = 0;
    Max_OneD_Grids = 0;
    //if (2<=level_stdout) printf("\n***** UCell_Box(MD_iter,1,CpyCell) *****\n");

    UCell_Box(1,1,CpyCell);

    My_Max_GridN_Atom = Max_GridN_Atom;

    MPI_Reduce(&My_Max_GridN_Atom, &Max_GridN_Atom, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_GridN_Atom, 1, MPI_INT, Host_ID, mpi_comm_level1);

    My_Max_OneD_Grids = Max_OneD_Grids;
    MPI_Reduce(&My_Max_OneD_Grids, &Max_OneD_Grids, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_OneD_Grids, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);

    //if (2<=level_stdout){
    //  printf("Max_OneD_Grids=%2d\n",Max_OneD_Grids);
    //}

    /*************************************
      find 
            Max_NumOLG
    *************************************/

    Max_NumOLG = 0;
    //if (2<=level_stdout) printf("\n***** UCell_Box(MD_iter,2,CpyCell) *****\n");

    UCell_Box(1,2,CpyCell);

    My_Max_NumOLG = Max_NumOLG;
    MPI_Reduce(&My_Max_NumOLG, &Max_NumOLG, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_NumOLG, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[12] = (int)(Max_NumOLG*ScaleSize) + 1;

    //if (2<=level_stdout){
    //  printf("YOUSO11=%i YOUSO12=%i YOUSO17=%i\n",
    //         List_YOUSO[11],List_YOUSO[12],List_YOUSO[17]);
    //}

    /*************************************
      setting of
                 GListTAtoms1
                 GListTAtoms2
    *************************************/

    //if (myid==Host_ID && 0<level_stdout)  printf("<UCell_Box> Info. of cutoff energy and num. of grids\n");

    UCell_Box(1,0,CpyCell);

    My_Max_GridN_Atom = Max_GridN_Atom;
    MPI_Reduce(&My_Max_GridN_Atom, &Max_GridN_Atom, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_GridN_Atom, 1, MPI_INT, Host_ID, mpi_comm_level1);

    My_Max_OneD_Grids = Max_OneD_Grids;
    MPI_Reduce(&My_Max_OneD_Grids, &Max_OneD_Grids, 1,
                MPI_INT, MPI_MAX, Host_ID, mpi_comm_level1);
    MPI_Bcast(&Max_OneD_Grids, 1, MPI_INT, Host_ID, mpi_comm_level1);
    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);
//  }
//
//  else if (UCell_flag==1) {
//
//    if (myid==Host_ID && 0<level_stdout) printf("<UCell_Box> Info. of cutoff energy and num. of grids\n");
//
//    UCell_Box(MD_iter,0,CpyCell);
//
//    List_YOUSO[11] = (int)(Max_GridN_Atom*ScaleSize) + 1;
//    List_YOUSO[17] = (int)(Max_OneD_Grids*ScaleSize);
//  }

  /****************************************************
     allocation of arrays:
  ****************************************************/

    /* arrays for the partitions B and C */

    if (SpinP_switch==3){ /* spin non-collinear */
      Density_Grid_B = (double**)malloc(sizeof(double*)*4); 
      for (k=0; k<=3; k++){
        Density_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
        for (i=0; i<My_NumGridB_AB; i++) Density_Grid_B[k][i] = 0.0;
      }
      /* added by fukuda 30.Dec.2014 */
      dDensity_Grid_B = (double***)malloc(sizeof(double**)*3); 
      for (k=0; k<3; k++){
         dDensity_Grid_B[k] = (double**)malloc(sizeof(double*)*6); 
         for (j=0; j<=5; j++){
           dDensity_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dDensity_Grid_B[k][j][i] = 0.0;
         }
      }
      /* fukuda */
      /* added by fukuda 24.Jan.2015 */
      ddDensity_Grid_B = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         ddDensity_Grid_B[k] = (double**)malloc(sizeof(double*)*6); 
         for (j=0; j<=5; j++){
           ddDensity_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) ddDensity_Grid_B[k][j][i] = 0.0;
         }
      }
      /* added by fukuda 2.May.2015 */
      dd2Density_Grid_B = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         dd2Density_Grid_B[k] = (double**)malloc(sizeof(double*)*6); 
         for (j=0; j<=5; j++){
           dd2Density_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dd2Density_Grid_B[k][j][i] = 0.0;
         }
      }
      /* fukuda */
    }
    else{
      Density_Grid_B = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Density_Grid_B[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
        for (i=0; i<My_NumGridB_AB; i++) Density_Grid_B[k][i] = 0.0;
      }
      /* added by fukuda 24.Jan.2015 */
      dDensity_Grid_B = (double***)malloc(sizeof(double**)*3); 
      for (k=0; k<3; k++){
         dDensity_Grid_B[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           dDensity_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dDensity_Grid_B[k][j][i] = 0.0;
         }
      }
      /* fukuda */
      /* added by fukuda 24.Jan.2015 */
      ddDensity_Grid_B = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         ddDensity_Grid_B[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           ddDensity_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) ddDensity_Grid_B[k][j][i] = 0.0;
         }
      }
      /* added by fukuda 2.May.2015 */
      dd2Density_Grid_B = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         dd2Density_Grid_B[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           dd2Density_Grid_B[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dd2Density_Grid_B[k][j][i] = 0.0;
         }
      }
      /* fukuda */
    }

      /* added by fukuda 23.Jun.2015 */
      /*** for NEGF ***/
      Density_Grid_B_i = (double**)malloc(sizeof(double*)*2); 
      for (k=0; k<=1; k++){
        Density_Grid_B_i[k] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
        for (i=0; i<My_NumGridB_AB; i++) Density_Grid_B_i[k][i] = 0.0;
      }
      dDensity_Grid_B_i = (double***)malloc(sizeof(double**)*3); 
      for (k=0; k<3; k++){
         dDensity_Grid_B_i[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           dDensity_Grid_B_i[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dDensity_Grid_B_i[k][j][i] = 0.0;
         }
      }
      ddDensity_Grid_B_i = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         ddDensity_Grid_B_i[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           ddDensity_Grid_B_i[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) ddDensity_Grid_B_i[k][j][i] = 0.0;
         }
      }
      dd2Density_Grid_B_i = (double***)malloc(sizeof(double**)*7); 
      for (k=0; k<7; k++){
         dd2Density_Grid_B_i[k] = (double**)malloc(sizeof(double*)*2); 
         for (j=0; j<=1; j++){
           dd2Density_Grid_B_i[k][j] = (double*)malloc(sizeof(double)*My_NumGridB_AB); 
           for (i=0; i<My_NumGridB_AB; i++) dd2Density_Grid_B_i[k][j][i] = 0.0;
         }
      }
      /* end for NEGF fukuda */

    /* Orbs_Grid */
    size_Orbs_Grid = 0;
    Orbs_Grid = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+1)); 
    Orbs_Grid[0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    Orbs_Grid[0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = F_M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      /* AITUNE */
      Orbs_Grid[Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*GridN_Atom[Gc_AN]); 
      int Nc;
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
        Orbs_Grid[Mc_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Spe_Total_NO[Cwan]); 
	size_Orbs_Grid += Spe_Total_NO[Cwan];
      }
      /* AITUNE */
    }

    /* dOrbs_Grid added by fukuda 30.Dec.2014 */
    size_dOrbs_Grid = 0;
    dOrbs_Grid = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*3); 
    for (k=0; k<3; k++){
      dOrbs_Grid[k] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+1)); 
      dOrbs_Grid[k][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
      dOrbs_Grid[k][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        /* AITUNE */
        dOrbs_Grid[k][Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*GridN_Atom[Gc_AN]); 
        int Nc;
        for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
          dOrbs_Grid[k][Mc_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Spe_Total_NO[Cwan]); 
	       size_dOrbs_Grid += Spe_Total_NO[Cwan];
        }
        /* AITUNE */
      }
    }

    //int Gc_AN2, Cwan2, NO0, NO02;
    //for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    //  Gc_AN = F_M2G[Mc_AN];
    //  Gc_AN2 = M2G[Mc_AN];
    //  Cwan = WhatSpecies[Gc_AN];
    //  Cwan2 = WhatSpecies[Gc_AN2];
    //  NO0 = Spe_Total_CNO[Cwan]; 
    //  NO02 = Spe_Total_CNO[Cwan2]; 
    //  printf("Mc_AN, Gc_AN, Gc_AN2, NO0, NO02, GridN_Atom, GridN_Atom2, = %d, %d, %d, %d, %d, %d, %d\n",Mc_AN, Gc_AN, Gc_AN2, NO0, NO02, GridN_Atom[Gc_AN], GridN_Atom[Gc_AN2] );fflush(stdout);
    //}

    /* ddOrbs_Grid added by fukuda 24.Jan.2015 */
    size_ddOrbs_Grid = 0;
    ddOrbs_Grid = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*7); 
    for (k=0; k<7; k++){
      ddOrbs_Grid[k] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+1)); 
      ddOrbs_Grid[k][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
      ddOrbs_Grid[k][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        /* AITUNE */
        ddOrbs_Grid[k][Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*GridN_Atom[Gc_AN]); 
        int Nc;
        for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
          ddOrbs_Grid[k][Mc_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*Spe_Total_NO[Cwan]); 
	       size_ddOrbs_Grid += Spe_Total_NO[Cwan];
        }
        /* AITUNE */
      }
    }

    /* COrbs_Grid */
    size_COrbs_Grid = 0;
    if (Cnt_switch!=0){
      COrbs_Grid = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+MatomnumF+1)); 
      COrbs_Grid[0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
      COrbs_Grid[0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
      for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
        Gc_AN = F_M2G[Mc_AN];
        Cwan = WhatSpecies[Gc_AN];
        COrbs_Grid[Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*Spe_Total_CNO[Cwan]); 
        for (i=0; i<Spe_Total_CNO[Cwan]; i++){
          COrbs_Grid[Mc_AN][i] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*GridN_Atom[Gc_AN]); 
          size_COrbs_Grid += GridN_Atom[Gc_AN];
        }
      }
    }

    /* CdOrbs_Grid added by fukuda 30.Dec.2014 */
    size_CdOrbs_Grid = 0;
    if (Cnt_switch!=0){
      CdOrbs_Grid = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*3); 
      for (k=0; k<3; k++){
        CdOrbs_Grid[k] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+MatomnumF+1)); 
        CdOrbs_Grid[k][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
        CdOrbs_Grid[k][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
        for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          CdOrbs_Grid[k][Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*Spe_Total_CNO[Cwan]); 
          for (i=0; i<Spe_Total_CNO[Cwan]; i++){
            CdOrbs_Grid[k][Mc_AN][i] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*GridN_Atom[Gc_AN]); 
            size_CdOrbs_Grid += GridN_Atom[Gc_AN];
          }
        }
      }
    }

    /* CddOrbs_Grid added by fukuda 30.Dec.2014 */
    size_CddOrbs_Grid = 0;
    if (Cnt_switch!=0){
      CddOrbs_Grid = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*7); 
      for (k=0; k<7; k++){
        CddOrbs_Grid[k] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(Matomnum+MatomnumF+1)); 
        CddOrbs_Grid[k][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
        CddOrbs_Grid[k][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
        for (Mc_AN=1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
          Gc_AN = F_M2G[Mc_AN];
          Cwan = WhatSpecies[Gc_AN];
          CddOrbs_Grid[k][Mc_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*Spe_Total_CNO[Cwan]); 
          for (i=0; i<Spe_Total_CNO[Cwan]; i++){
            CddOrbs_Grid[k][Mc_AN][i] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*GridN_Atom[Gc_AN]); 
            size_CddOrbs_Grid += GridN_Atom[Gc_AN];
          }
        }
      }
    }

    /* Orbs_Grid_FNAN */
    size_Orbs_Grid_FNAN = 0;
    Orbs_Grid_FNAN = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*(Matomnum+1)); 
    Orbs_Grid_FNAN[0] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*1); 
    Orbs_Grid_FNAN[0][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    Orbs_Grid_FNAN[0][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];    
      Orbs_Grid_FNAN[Mc_AN] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(FNAN[Gc_AN]+1)); 

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

	  Mh_AN = F_G2M[Gh_AN];
	  Rnh = ncn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          NO1 = Spe_Total_NO[Hwan];
          /* AITUNE */
          if (0<NumOLG[Mc_AN][h_AN]){
  	    Orbs_Grid_FNAN[Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*NumOLG[Mc_AN][h_AN]); 
            int Nc;
	    for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
              Orbs_Grid_FNAN[Mc_AN][h_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*NO1); 
              size_Orbs_Grid_FNAN += NO1;
	    }
	  }
          /* AITUNE */
	}
        else {
          Orbs_Grid_FNAN[Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
          Orbs_Grid_FNAN[Mc_AN][h_AN][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
          size_Orbs_Grid_FNAN += 1;
        }
      }
    }

    /* dOrbs_Grid_FNAN added by fukuda 30.Dec.2014 */
    size_dOrbs_Grid_FNAN = 0;
    dOrbs_Grid_FNAN = (Type_Orbs_Grid*****)malloc(sizeof(Type_Orbs_Grid****)*3); 
    for (k=0; k<3; k++){
    dOrbs_Grid_FNAN[k] = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*(Matomnum+1)); 
    dOrbs_Grid_FNAN[k][0] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*1); 
    dOrbs_Grid_FNAN[k][0][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    dOrbs_Grid_FNAN[k][0][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];    
      dOrbs_Grid_FNAN[k][Mc_AN] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(FNAN[Gc_AN]+1)); 

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

	  Mh_AN = F_G2M[Gh_AN];
	  Rnh = ncn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          NO1 = Spe_Total_NO[Hwan];
          /* AITUNE */
          if (0<NumOLG[Mc_AN][h_AN]){
  	    dOrbs_Grid_FNAN[k][Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*NumOLG[Mc_AN][h_AN]); 
            int Nc;
	    for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
              dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*NO1); 
              size_dOrbs_Grid_FNAN += NO1;
	    }
	  }
          /* AITUNE */
	}
        else {
          dOrbs_Grid_FNAN[k][Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
          dOrbs_Grid_FNAN[k][Mc_AN][h_AN][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
          size_dOrbs_Grid_FNAN += 1;
        }
      }
    }
    } /* k */

    /* ddOrbs_Grid_FNAN added by fukuda 30.Dec.2014 */
    size_ddOrbs_Grid_FNAN = 0;
    ddOrbs_Grid_FNAN = (Type_Orbs_Grid*****)malloc(sizeof(Type_Orbs_Grid****)*7); 
    for (k=0; k<7; k++){
    ddOrbs_Grid_FNAN[k] = (Type_Orbs_Grid****)malloc(sizeof(Type_Orbs_Grid***)*(Matomnum+1)); 
    ddOrbs_Grid_FNAN[k][0] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*1); 
    ddOrbs_Grid_FNAN[k][0][0] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
    ddOrbs_Grid_FNAN[k][0][0][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = M2G[Mc_AN];    
      ddOrbs_Grid_FNAN[k][Mc_AN] = (Type_Orbs_Grid***)malloc(sizeof(Type_Orbs_Grid**)*(FNAN[Gc_AN]+1)); 

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

        Gh_AN = natn[Gc_AN][h_AN];

        if (G2ID[Gh_AN]!=myid){

	  Mh_AN = F_G2M[Gh_AN];
	  Rnh = ncn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
          NO1 = Spe_Total_NO[Hwan];
          /* AITUNE */
          if (0<NumOLG[Mc_AN][h_AN]){
  	    ddOrbs_Grid_FNAN[k][Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*NumOLG[Mc_AN][h_AN]); 
            int Nc;
	    for (Nc=0; Nc<NumOLG[Mc_AN][h_AN]; Nc++){
              ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nc] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*NO1); 
              size_ddOrbs_Grid_FNAN += NO1;
	    }
	  }
          /* AITUNE */
	}
        else {
          ddOrbs_Grid_FNAN[k][Mc_AN][h_AN] = (Type_Orbs_Grid**)malloc(sizeof(Type_Orbs_Grid*)*1); 
          ddOrbs_Grid_FNAN[k][Mc_AN][h_AN][0] = (Type_Orbs_Grid*)malloc(sizeof(Type_Orbs_Grid)*1); 
          size_ddOrbs_Grid_FNAN += 1;
        }
      }
    }
    } /* k */
  
}






void UCell_Box(int MD_iter, int estimate_switch, int CpyCell)
{
  static int firsttime=1;
  int size_GListTAtoms1;
  int po,N3[4],i,j,k;
  int size_GridListAtom,size_MGridListAtom;
  int NOC[4],nn1,nn2,nn3,l1,l2,l3,lmax;
  int ct_AN,N,n1,n2,n3,Cwan,Rn,Rn1,Ng1,Ng2,Ng3;
  int p,q,r,s,pmax,qmax,rmax,smax,popt,qopt,ropt,sopt;
  int Nct,MinN,Scale,ScaleA,ScaleB,ScaleC,Nm[4];
  int Mc_AN,Mc_AN0,Gc_AN,h_AN,h_AN0,Mh_AN,Gh_AN,Nh,Rh,Nog,GNh,GRh;
  int ll1,ll2,ll3,Nnb,My_Max_NumOLG;
  int lll1,lll2,lll3,GRh1,Nc,GNc,GRc,size_array;
  int *TAtoms0,*TCells0,*TAtoms1,*TAtoms2;
  int **Tmp_GridListAtom,**Tmp_CellListAtom;
  int nmin[4],nmax[4],Np;
  double LgTN,LgN,Lg2,Lg3,Lg5,Lg7,DouN[4],A1,A2,A3;
  double MinD,MinR,CutR2,r2,coef;
  double sa,sa_cri,tmp[4],Cxyz[4];
  double b[4],c[4],v[4],rcut;
  double xc,yc,zc,xm,ym,zm;
  double dx,dy,dz,sn1,sn2,sn3;
  double B2,C2,CellV;
  double S_Lng,L_Lng,LngA,LngB,LngC,x,y,z;
  double GVolume,buffer_scale,GridV;
  double stime,etime;

  int *TempGrid,*TempCell;
  int numprocs,myid,ID,IDS,IDR,tag=999;
  double Stime_atom, Etime_atom;
  char file_UCell[YOUSO10];
  FILE *fp;
  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;
  int alloc_first0, alloc_first2;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);



  /****************************************************
                Reciprocal lattice vectors
  ****************************************************/

  if (estimate_switch<=1){
  
    Cross_Product(tv[2],tv[3],tmp);
    CellV = Dot_Product(tv[1],tmp); 
  
    Cross_Product(tv[2],tv[3],tmp);
    rtv[1][1] = 2.0*PI*tmp[1]/CellV;
    rtv[1][2] = 2.0*PI*tmp[2]/CellV;
    rtv[1][3] = 2.0*PI*tmp[3]/CellV;
  
    Cross_Product(tv[3],tv[1],tmp);
    rtv[2][1] = 2.0*PI*tmp[1]/CellV;
    rtv[2][2] = 2.0*PI*tmp[2]/CellV;
    rtv[2][3] = 2.0*PI*tmp[3]/CellV;
  
    Cross_Product(tv[1],tv[2],tmp);
    rtv[3][1] = 2.0*PI*tmp[1]/CellV;
    rtv[3][2] = 2.0*PI*tmp[2]/CellV;
    rtv[3][3] = 2.0*PI*tmp[3]/CellV;

    //if (myid==Host_ID && 0<level_stdout){    

    //  printf("lattice vectors (bohr)\n");
    //  printf("A  = %15.12f, %15.12f, %15.12f\n",tv[1][1],tv[1][2],tv[1][3]);
    //  printf("B  = %15.12f, %15.12f, %15.12f\n",tv[2][1],tv[2][2],tv[2][3]);
    //  printf("C  = %15.12f, %15.12f, %15.12f\n",tv[3][1],tv[3][2],tv[3][3]);

    //  printf("reciprocal lattice vectors (bohr^-1)\n");
    //  printf("RA = %15.12f, %15.12f, %15.12f\n",rtv[1][1],rtv[1][2],rtv[1][3]);
    //  printf("RB = %15.12f, %15.12f, %15.12f\n",rtv[2][1],rtv[2][2],rtv[2][3]);
    //  printf("RB = %15.12f, %15.12f, %15.12f\n",rtv[3][1],rtv[3][2],rtv[3][3]);
    //}  

    /* calculate gtv, rgtv, A2, B2, C2, and used cutoff energies */

    gtv[1][1] = tv[1][1]/(double)Ngrid1;
    gtv[1][2] = tv[1][2]/(double)Ngrid1;
    gtv[1][3] = tv[1][3]/(double)Ngrid1;
    
    gtv[2][1] = tv[2][1]/(double)Ngrid2;
    gtv[2][2] = tv[2][2]/(double)Ngrid2;
    gtv[2][3] = tv[2][3]/(double)Ngrid2;
    
    gtv[3][1] = tv[3][1]/(double)Ngrid3;
    gtv[3][2] = tv[3][2]/(double)Ngrid3;
    gtv[3][3] = tv[3][3]/(double)Ngrid3;

    Cross_Product(gtv[2],gtv[3],tmp);

    GridV = Dot_Product(gtv[1],tmp); 
    GVolume = fabs( GridV );

    Cross_Product(gtv[2],gtv[3],tmp);
    rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
    rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
    rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

    Cross_Product(gtv[3],gtv[1],tmp);
    rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
    rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
    rgtv[2][3] = 2.0*PI*tmp[3]/GridV;
    
    Cross_Product(gtv[1],gtv[2],tmp);
    rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
    rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
    rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

    A2 = rgtv[1][1]*rgtv[1][1] + rgtv[1][2]*rgtv[1][2] + rgtv[1][3]*rgtv[1][3];
    B2 = rgtv[2][1]*rgtv[2][1] + rgtv[2][2]*rgtv[2][2] + rgtv[2][3]*rgtv[2][3];
    C2 = rgtv[3][1]*rgtv[3][1] + rgtv[3][2]*rgtv[3][2] + rgtv[3][3]*rgtv[3][3];

    A2 = A2/4.0;  /* note: change the unit from Hatree to Rydberg by multiplying 1/2 */
    B2 = B2/4.0;
    C2 = C2/4.0;

    /* print information to std output */

    //if (estimate_switch==0 || 2<=level_stdout){
    //  if (myid==Host_ID && 0<level_stdout) {
    //    printf("    Used cutoff energy (Ryd) for 3D-grids = %7.4f, %7.4f, %7.4f\n",
    //            A2,B2,C2);
    //    printf("Num. of grids of a-, b-, and c-axes = %2d, %2d, %2d\n",
    //            Ngrid1,Ngrid2,Ngrid3);
    //  }
    //}

    Max_OneD_Grids = Ngrid1;
    if (Max_OneD_Grids<Ngrid2) Max_OneD_Grids = Ngrid2;
    if (Max_OneD_Grids<Ngrid3) Max_OneD_Grids = Ngrid3;

  } /* if (estimate_switch<=1) */

  /****************************************************
       Setting the center of unit cell and grids
  ****************************************************/

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

  /* added by T.Ohwaki */
  if(length_gtv[1] != 0.0 && ESM_switch!=0){

    double xc_tmp; /* added by T.Ohwaki */

    xc_tmp = (int)(xc / length_gtv[1]);
    xc = ((double)xc_tmp * length_gtv[1]);
    //if(myid==Host_ID && 0<level_stdout){
    //  printf("<ESM> length_gtv[1],xc/length_gtv[1] = %12.9f,%12.9f \n",length_gtv[1],xc/length_gtv[1]);
    //}
  } /* added by T.Ohwaki */

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

    //if ( 1.0e+8<scf_fixed_dens_origin[0] &&
    //     1.0e+8<scf_fixed_dens_origin[1] &&
	 //1.0e+8<scf_fixed_dens_origin[2] ){

    //  Grid_Origin[1] = xc - xm;
    //  Grid_Origin[2] = yc - ym;
    //  Grid_Origin[3] = zc - zm;
    //  ///*** fukuda temp ***/
    //  //Grid_Origin[1] = 31.18048; /* Graphene arm-chair*/
    //  //Grid_Origin[1] = 13.83279; /* Graphene ZGNR HCCCCH*/
    //  //Grid_Origin[1] = 16.13826; /* Graphene ZGNR CCCC*/
    //  ///*** ***/

    //  /* added by T.Ohwaki */
    //  if (myid==Host_ID && ESM_switch!=0 && 0<level_stdout){
    //    printf("xc=%15.12f yc=%15.12f zc=%15.12f\n",xc,yc,zc);
    //    printf("xm=%15.12f ym=%15.12f zm=%15.12f\n",xm,ym,zm);
    //  }

    //}
    //else{
    //     Grid_Origin[1] = scf_fixed_dens_origin[0];
    //     Grid_Origin[2] = scf_fixed_dens_origin[1];
    //     Grid_Origin[3] = scf_fixed_dens_origin[2];
    //}
    
    //if ( 1.0e+8<scf_fixed_dens_origin[0] ){
    //  Grid_Origin[1] = xc - xm;
    //}
    //else{
    //     Grid_Origin[1] = scf_fixed_dens_origin[0];
    //}
    //if ( 1.0e+8<scf_fixed_dens_origin[1] ){
    //  Grid_Origin[2] = yc - ym;
    //}
    //else{
    //     Grid_Origin[2] = scf_fixed_dens_origin[1];
    //}
    //if ( 1.0e+8<scf_fixed_dens_origin[2] ){
    //  Grid_Origin[3] = zc - zm;
    //}
    //else{
    //     Grid_Origin[3] = scf_fixed_dens_origin[2];
    //}

      /* added by T.Ohwaki */
      //if (myid==Host_ID && ESM_switch!=0 && 0<level_stdout){
      //  printf("xc=%15.12f yc=%15.12f zc=%15.12f\n",xc,yc,zc);
      //  printf("xm=%15.12f ym=%15.12f zm=%15.12f\n",xm,ym,zm);
      //}
    
    //if (myid==Host_ID && 0<level_stdout){
    //  //printf("scf_fixed_dens_origin %15.12f %15.12f %15.12f\n",
    //  //        scf_fixed_dens_origin[0],scf_fixed_dens_origin[1],scf_fixed_dens_origin[2]);
    //  //printf("xc yc zc %15.12f %15.12f %15.12f\n",
    //  //        xc, yc, zc);
    //  //printf("xm ym zm %15.12f %15.12f %15.12f\n",
    //  //        xm, ym, zm);
    //  printf("Grid_Origin %15.12f %15.12f %15.12f\n",
    //          Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
    //}    

    //TNumGrid = Ngrid1*Ngrid2*Ngrid3;
    //Last_TNumGrid = (int)(ScaleSize*Ngrid1*Ngrid2*Ngrid3);

    //LastBoxCenterX = xc;
    //LastBoxCenterY = yc;
    //LastBoxCenterZ = zc;

  /****************************************************
            xyz-coordinate to cell-coordinate
  ****************************************************/

  //if (estimate_switch<=1){

  //  #pragma omp parallel shared(level_stdout,rtv,Cell_Gxyz,Grid_Origin,Gxyz,M2G,Matomnum) private(Mc_AN,Gc_AN,Cxyz,OMPID,Nthrds,Nprocs)
  //  {

  //    /* get info. on OpenMP */ 
  //
  //    OMPID = omp_get_thread_num();
  //    Nthrds = omp_get_num_threads();
  //    Nprocs = omp_get_num_procs();

  //    for (Mc_AN=(OMPID+1); Mc_AN<=Matomnum; Mc_AN+=Nthrds){

  //      Gc_AN = M2G[Mc_AN];
  //      Cxyz[1] = Gxyz[Gc_AN][1] - Grid_Origin[1];
  //      Cxyz[2] = Gxyz[Gc_AN][2] - Grid_Origin[2];
  //      Cxyz[3] = Gxyz[Gc_AN][3] - Grid_Origin[3];
  //      Cell_Gxyz[Gc_AN][1] = Dot_Product(Cxyz,rtv[1])*0.5/PI;
  //      Cell_Gxyz[Gc_AN][2] = Dot_Product(Cxyz,rtv[2])*0.5/PI;
  //      Cell_Gxyz[Gc_AN][3] = Dot_Product(Cxyz,rtv[3])*0.5/PI;

  //      if (2<=level_stdout){
  //        printf("Cell_Gxyz %3d  %15.12f %15.12f %15.12f\n",
  //      	 Gc_AN,
  //      	 Cell_Gxyz[Gc_AN][1],
  //      	 Cell_Gxyz[Gc_AN][2],
  //      	 Cell_Gxyz[Gc_AN][3]); 
  //      }
  //    }

  //  } /* #pragma omp parallel */

  //  /****************
  //   MPI: 
  //       Cell_Gxyz
  //  *****************/

  //  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
  //    ID = G2ID[ct_AN];
  //    MPI_Bcast(&Cell_Gxyz[ct_AN][0], 4, MPI_DOUBLE, ID, mpi_comm_level1);
  //  }
  //}

  /****************************************************
            Find grids overlaping to each atom
  ****************************************************/

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

  Cross_Product(gtv[2],gtv[3],tmp);
  GridV = Dot_Product(gtv[1],tmp);

  length_gtv[1] = sqrt( Dot_Product(gtv[1], gtv[1]) );
  length_gtv[2] = sqrt( Dot_Product(gtv[2], gtv[2]) );
  length_gtv[3] = sqrt( Dot_Product(gtv[3], gtv[3]) );

  /* rgtv */

  Cross_Product(gtv[2],gtv[3],tmp);
  rgtv[1][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[1][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[1][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv[3],gtv[1],tmp);
  rgtv[2][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[2][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[2][3] = 2.0*PI*tmp[3]/GridV;

  Cross_Product(gtv[1],gtv[2],tmp);
  rgtv[3][1] = 2.0*PI*tmp[1]/GridV;
  rgtv[3][2] = 2.0*PI*tmp[2]/GridV;
  rgtv[3][3] = 2.0*PI*tmp[3]/GridV;

  if ( (estimate_switch==0 || 2<=level_stdout) && myid==Host_ID ){

    if (0<level_stdout){ 

      printf("Cell vectors (bohr) of the grid cell (gtv)\n");
      printf("  gtv_a = %15.12f, %15.12f, %15.12f\n",gtv[1][1],gtv[1][2],gtv[1][3]);
      printf("  gtv_b = %15.12f, %15.12f, %15.12f\n",gtv[2][1],gtv[2][2],gtv[2][3]);
      printf("  gtv_c = %15.12f, %15.12f, %15.12f\n",gtv[3][1],gtv[3][2],gtv[3][3]);
      printf("  |gtv_a| = %15.12f\n",length_gtv[1]);
      printf("  |gtv_b| = %15.12f\n",length_gtv[2]);
      printf("  |gtv_c| = %15.12f\n",length_gtv[3]);
    }
  }

  /**********************************
    allocation of arrays: 

    Tmp_GridListAtom
    Tmp_CellListAtom
    MGridListAtom
  **********************************/
 
  Tmp_GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  Tmp_CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+MatomnumF+1));
  MGridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  Tmp_GridListAtom[0] = (int*)malloc(sizeof(int)*1);
  Tmp_CellListAtom[0] = (int*)malloc(sizeof(int)*1);
  MGridListAtom[0] = (int*)malloc(sizeof(int)*1);
  alloc_first2 = 0;

  /****************************************************
   1) find a neighbouring point of the atom Mc_AN
   2) the ranges which deterinie a box on the atom Mc_AN 
   3) determine whether overlap exists or not
  ****************************************************/

  /* for allocation of arrays */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    rcut = Spe_Atom_Cut1[Cwan] + 0.5;

    for (k=1; k<=3; k++){

      if      (k==1){ i = 2; j = 3; }
      else if (k==2){ i = 3; j = 1; }
      else if (k==3){ i = 1; j = 2; }

      b[1] = tv_ori[i][1];
      b[2] = tv_ori[i][2];
      b[3] = tv_ori[i][3];

      c[1] = tv_ori[j][1];
      c[2] = tv_ori[j][2];
      c[3] = tv_ori[j][3];

      Cross_Product(b,c,v);
      coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));

      v[1] = coef*v[1];
      v[2] = coef*v[2];
      v[3] = coef*v[3];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      /* find the maximum range of grids */
      nmax[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      /* find the mimum range of grids */
      nmin[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      if (nmax[k]<nmin[k]){
        i = nmin[k];
        j = nmax[k];
        nmin[k] = j;
        nmax[k] = i;
      } 
  
    } /* k */  

    /* allocation of arrays */ 

    Np = (nmax[1]-nmin[1]+1)*(nmax[2]-nmin[2]+1)*(nmax[3]-nmin[3]+1)*3/2;
    
    Tmp_GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);
    Tmp_CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);
    MGridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*Np);

  } /* Mc_AN */

  /* store Tmp_GridListAtom and Tmp_CellListAtom */

  size_array = (int)(Max_GridN_Atom*ScaleSize);
  TempGrid = (int*)malloc(sizeof(int)*size_array); 
  TempCell = (int*)malloc(sizeof(int)*size_array); 

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    Cwan = WhatSpecies[Gc_AN];
    rcut = Spe_Atom_Cut1[Cwan] + 0.5;

    for (k=1; k<=3; k++){

      if      (k==1){ i = 2; j = 3; }
      else if (k==2){ i = 3; j = 1; }
      else if (k==3){ i = 1; j = 2; }

      b[1] = tv_ori[i][1];
      b[2] = tv_ori[i][2];
      b[3] = tv_ori[i][3];

      c[1] = tv_ori[j][1];
      c[2] = tv_ori[j][2];
      c[3] = tv_ori[j][3];

      Cross_Product(b,c,v);
      coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));

      v[1] = coef*v[1];
      v[2] = coef*v[2];
      v[3] = coef*v[3];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      /* find the maximum range of grids */
      nmax[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      /* find the mimum range of grids */
      nmin[k] = Dot_Product(Cxyz,rgtv[k])*0.5/PI;
  
      if (nmax[k]<nmin[k]){
        i = nmin[k];
        j = nmax[k];
        nmin[k] = j;
        nmax[k] = i;
      } 

    } /* k */  

    CutR2 = Spe_Atom_Cut1[Cwan]*Spe_Atom_Cut1[Cwan];

    Nct = 0;
    for (n1=nmin[1]; n1<=nmax[1]; n1++){
      for (n2=nmin[2]; n2<=nmax[2]; n2++){
	for (n3=nmin[3]; n3<=nmax[3]; n3++){

	  Find_CGrids(1,n1,n2,n3,Cxyz,NOC);
	  Rn = NOC[0];
	  l1 = NOC[1];
	  l2 = NOC[2];
	  l3 = NOC[3];
	  N = l1*Ngrid2*Ngrid3 + l2*Ngrid3 + l3;

	  dx = Cxyz[1] - Gxyz[Gc_AN][1];
	  dy = Cxyz[2] - Gxyz[Gc_AN][2];
	  dz = Cxyz[3] - Gxyz[Gc_AN][3];

	  r2 = dx*dx + dy*dy + dz*dz;
	  if (r2<=CutR2){
	    if (estimate_switch!=1){
	      TempGrid[Nct+1] = N;
	      TempCell[Nct+1] = Rn;
	    }             
	    Nct++;
	  }
	}
      }
    }

    Np = (nmax[1]-nmin[1]+1)*(nmax[2]-nmin[2]+1)*(nmax[3]-nmin[3]+1)*3/2;
    if (Np<Nct){
      printf("Invalid access in truncation.c\n"); 
      MPI_Finalize();
      exit(0); 
    }

    GridN_Atom[Gc_AN] = Nct;

    if (estimate_switch!=1){

      /* sorting */
      qsort_int((long)Nct,TempGrid,TempCell);

      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	Tmp_GridListAtom[Mc_AN][Nc] = TempGrid[Nc+1];
	Tmp_CellListAtom[Mc_AN][Nc] = TempCell[Nc+1];
      }
    }
  }

  free(TempCell);
  free(TempGrid);

  /* calculate size_GridListAtom */

  size_GridListAtom = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    size_GridListAtom += GridN_Atom[Gc_AN];
    if (Max_GridN_Atom<GridN_Atom[Gc_AN]) Max_GridN_Atom = GridN_Atom[Gc_AN];
  }

  /****************************************************
   MPI: 

       GridN_Atom
  ****************************************************/

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    ID = G2ID[ct_AN];
    MPI_Bcast(&GridN_Atom[ct_AN], 1, MPI_INT, ID, mpi_comm_level1);
  }

  //if (myid==Host_ID && estimate_switch==0 && 0<level_stdout){
  //  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
  //    if (ct_AN<=20 && level_stdout<=1){
  //       printf("Num. of grids overlapping with atom %4d = %4d\n",
  //               ct_AN, GridN_Atom[ct_AN]);
  //    }
  //  }

  //  if (20<atomnum && level_stdout<=1){
  //    printf("     ..........\n");
  //    printf("     ......\n\n");
  //  }
  //}

  /****************************************************
    allocation of arrays:

    Tmp_GridListAtom
    Tmp_CellListAtom
  ****************************************************/
  
  size_MGridListAtom = size_GridListAtom;

  for (Mc_AN=Matomnum+1; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    Tmp_GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    Tmp_CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    size_MGridListAtom += GridN_Atom[Gc_AN];
  }

  /****************************************************
   MPI: 

       Tmp_GridListAtom
       Tmp_CellListAtom
  ****************************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  /* Tmp_GridListAtom */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      /* Sending of data to IDS */
      if (F_Snd_Num[IDS]!=0){
        for (i=0; i<F_Snd_Num[IDS]; i++){
          Mc_AN = Snd_MAN[IDS][i];
          Gc_AN = Snd_GAN[IDS][i];
          MPI_Isend(&Tmp_GridListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (F_Rcv_Num[IDR]!=0){
        tag = 999;
        for (i=0; i<F_Rcv_Num[IDR]; i++){
          Mc_AN = F_TopMAN[IDR] + i;
          Gc_AN = F_M2G[Mc_AN];
          MPI_Recv(&Tmp_GridListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	}          
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    }     
  }

  /* Tmp_CellListAtom */

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      /* Sending of data to IDS */
      if (F_Snd_Num[IDS]!=0){
        for (i=0; i<F_Snd_Num[IDS]; i++){
          Mc_AN = Snd_MAN[IDS][i];
          Gc_AN = Snd_GAN[IDS][i];
          MPI_Isend(&Tmp_CellListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                    IDS, tag, mpi_comm_level1, &request);
	}
      }

      /* Receiving of data from IDR */
      if (F_Rcv_Num[IDR]!=0){
        tag = 999;
        for (i=0; i<F_Rcv_Num[IDR]; i++){
          Mc_AN = F_TopMAN[IDR] + i;
          Gc_AN = F_M2G[Mc_AN];
          MPI_Recv(&Tmp_CellListAtom[Mc_AN][0], GridN_Atom[Gc_AN], MPI_INT,
                   IDR, tag, mpi_comm_level1, &stat);
	}          
      }

      if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
    }     
  }

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  if (estimate_switch!=1){

    /****************************************************
            Find overlap grids between two orbitals
    ****************************************************/
    
    if (estimate_switch==0){
      size_GListTAtoms1 = 0;

      GListTAtoms1 = (int***)malloc(sizeof(int**)*(Matomnum+1));
      GListTAtoms2 = (int***)malloc(sizeof(int**)*(Matomnum+1));
      alloc_first0 = 0;
    }

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];

      if (Mc_AN==0){

	FNAN[0] = 0; 

        if (estimate_switch==0){

          GListTAtoms1[0] = (int**)malloc(sizeof(int*)*1);
          GListTAtoms2[0] = (int**)malloc(sizeof(int*)*1);

          GListTAtoms1[0][0] = (int*)malloc(sizeof(int)*1);
          GListTAtoms2[0][0] = (int*)malloc(sizeof(int)*1);
        }
      }

      else{

        if (estimate_switch==0){

          GListTAtoms1[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
          GListTAtoms2[Mc_AN] = (int**)malloc(sizeof(int*)*(FNAN[Gc_AN]+1));
        }

        h_AN0 = 0;

#pragma omp parallel shared(List_YOUSO,GListTAtoms1,GListTAtoms2,ScaleSize,Max_NumOLG,size_GListTAtoms1,level_stdout,NumOLG,estimate_switch,CpyCell,Mc_AN,Tmp_CellListAtom,Tmp_GridListAtom,GridN_Atom,atv_ijk,ncn,F_G2M,natn,h_AN0,FNAN,Gc_AN) private(OMPID,Nthrds,Nprocs,h_AN,Gh_AN,Mh_AN,Rh,l1,l2,l3,Nog,Nc,Nh,GNh,GRh,ll1,ll2,ll3,lll1,lll2,lll3,GRh1,po,GNc,GRc,TAtoms0,TCells0,TAtoms1,TAtoms2,size_array,i)
	{        

	  /*******************************************************
            allocation of temporal arrays
	  *******************************************************/

	  size_array = (int)((Max_NumOLG+2)*ScaleSize);
	  TAtoms0 = (int*)malloc(sizeof(int)*size_array);
	  TCells0 = (int*)malloc(sizeof(int)*size_array);
	  TAtoms1 = (int*)malloc(sizeof(int)*size_array);
	  TAtoms2 = (int*)malloc(sizeof(int)*size_array);

	  /* get info. on OpenMP */ 
  
	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

          do {  

#pragma omp barrier
            h_AN = h_AN0 + OMPID;

            if (h_AN<=FNAN[Gc_AN]){

	      Gh_AN = natn[Gc_AN][h_AN];
	      Mh_AN = F_G2M[Gh_AN];
	      Rh = ncn[Gc_AN][h_AN];

	      l1 = atv_ijk[Rh][1];
	      l2 = atv_ijk[Rh][2];
	      l3 = atv_ijk[Rh][3];
  
	      Nog = -1;
	      Nc = 0;

	      for (Nh=0; Nh<GridN_Atom[Gh_AN]; Nh++){

		GNh = Tmp_GridListAtom[Mh_AN][Nh];
		GRh = Tmp_CellListAtom[Mh_AN][Nh];

		ll1 = atv_ijk[GRh][1];
		ll2 = atv_ijk[GRh][2];
		ll3 = atv_ijk[GRh][3];

		lll1 = l1 + ll1;
		lll2 = l2 + ll2;
		lll3 = l3 + ll3;

		if (Tmp_GridListAtom[Mc_AN][0]<=GNh) {

		  /* find the initial Nc */

		  if (GNh==0) {
		    Nc = 0;
		  }
		  else {
		    while ( GNh<=Tmp_GridListAtom[Mc_AN][Nc] && Nc!=0 ){
		      Nc = Nc - 10;
		      if (Nc<0) Nc = 0;
		    }
		  }

		  /*  find whether there is the overlapping or not. */

		  if (abs(lll1)<=CpyCell && abs(lll2)<=CpyCell && abs(lll3)<=CpyCell){

		    //GRh1 = R_atv(CpyCell,lll1,lll2,lll3);
                    GRh1 = ratv[lll1+CpyCell][lll2+CpyCell][lll3+CpyCell];

		    po = 0;

		    while (po==0 && Nc<GridN_Atom[Gc_AN]) {

		      GNc = Tmp_GridListAtom[Mc_AN][Nc];
		      GRc = Tmp_CellListAtom[Mc_AN][Nc];

		      if (GNc==GNh && GRc==GRh1){

			Nog++;

			if (estimate_switch==0){
			  TAtoms0[Nog] = GNc;
			  TCells0[Nog] = GRc;
			  TAtoms1[Nog] = Nc;
			  TAtoms2[Nog] = Nh;
			}

			po = 1;
		      }

		      else if (GNh<GNc){
			po = 1; 
		      } 

		      Nc++;

		    } /* while (...) */

		    /* for Nc==GridN_Atom[Gc_AN] */

		    Nc--;
		    if (Nc<0) Nc = 0;

		  } /* if (abs.... ) */
		} /* if (Tmp_GridListAtom[Mc_AN][0]<=GNh) */
	      } /* Nh */

	      NumOLG[Mc_AN][h_AN] = Nog + 1;

	      if (List_YOUSO[12]<(Nog+1) && estimate_switch==0){
		printf("YOUSO12<(Nog+1)\n");
                MPI_Finalize();
		exit(1);
	      }

	      if (2<=level_stdout){
		printf("Num. of grids overlapping between atoms %2d (G) and %2d (L) (%2d) = %2d\n",
		       Gc_AN,h_AN,Gh_AN,Nog+1);
	      }
	    
	    } /* if (h_AN<=FNAN[Gc_AN]) */

#pragma omp barrier
#pragma omp flush(NumOLG)

	    if (estimate_switch==0){

              /* allocation of arrays */

              if (OMPID==0){

		for (i=h_AN0; i<(h_AN0+Nthrds); i++){

                  if (i<=FNAN[Gc_AN]){

		    size_GListTAtoms1 += NumOLG[Mc_AN][i];

		    GListTAtoms1[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);
		    GListTAtoms2[Mc_AN][i] = (int*)malloc(sizeof(int)*NumOLG[Mc_AN][i]);

		  }
		}
	      }

#pragma omp barrier
#pragma omp flush(GListTAtoms1,GListTAtoms2)

              if (h_AN<=FNAN[Gc_AN]){

		for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){

		  GListTAtoms1[Mc_AN][h_AN][Nog] = TAtoms1[Nog];
		  GListTAtoms2[Mc_AN][h_AN][Nog] = TAtoms2[Nog];
		}
	      }
	    }

            /* increament of h_AN0 */

            if (OMPID==0) h_AN0 += Nthrds;
#pragma omp barrier
#pragma omp flush(h_AN0)

	  } while (h_AN0<=FNAN[Gc_AN]);

          /* freeing of arrays */

	  free(TAtoms2);
	  free(TAtoms1);
	  free(TCells0);
	  free(TAtoms0);

	} /* #pragma omp parallel */

      }

    } /* Mc_AN */ 

    /* find Max_NumOLG */

    My_Max_NumOLG = 0;
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
        if (My_Max_NumOLG<NumOLG[Mc_AN][h_AN]) My_Max_NumOLG = NumOLG[Mc_AN][h_AN];
      }
    }

    MPI_Allreduce(&My_Max_NumOLG, &Max_NumOLG, 1, MPI_INT, MPI_MAX, mpi_comm_level1);

  } /* if (estimate_switch!=1) */

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  /****************************************************
       Tmp_GridListAtom -> GridListAtom
       Tmp_CellListAtom -> CellListAtom                     
  ****************************************************/

  size_GridListAtom = 0;

  GridListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  GridListAtom[0] = (int*)malloc(sizeof(int)*1);
  size_GridListAtom++;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    GridListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    size_GridListAtom += GridN_Atom[Gc_AN];
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      GridListAtom[Mc_AN][Nc] = Tmp_GridListAtom[Mc_AN][Nc];
    }
  }

  CellListAtom = (int**)malloc(sizeof(int*)*(Matomnum+1));
  CellListAtom[0] = (int*)malloc(sizeof(int)*1);
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    CellListAtom[Mc_AN] = (int*)malloc(sizeof(int)*GridN_Atom[Gc_AN]);
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      CellListAtom[Mc_AN][Nc] = Tmp_CellListAtom[Mc_AN][Nc];
    }
  }

  /****************************************************
   construct the data structure for MPI communications
   for grid data
  ****************************************************/

  if (estimate_switch==0){

    Construct_MPI_Data_Structure_Grid();

    Ng1 = Max_Grid_Index[1] - Min_Grid_Index[1] + 1;
    Ng2 = Max_Grid_Index[2] - Min_Grid_Index[2] + 1;
    Ng3 = Max_Grid_Index[3] - Min_Grid_Index[3] + 1;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      Gc_AN = F_M2G[Mc_AN];

#pragma omp parallel shared(Min_Grid_Index,atv_ijk,Ng1,Ng2,Ng3,Ngrid1,Ngrid2,Ngrid3,MGridListAtom,Mc_AN,Tmp_GridListAtom,Tmp_CellListAtom,GridN_Atom,Gc_AN) private(Nc,GNc,GRc,N3,n1,n2,n3,N,OMPID,Nthrds,Nprocs)
      {

	/* get info. on OpenMP */ 
  
	OMPID = omp_get_thread_num();
	Nthrds = omp_get_num_threads();
	Nprocs = omp_get_num_procs();

	for (Nc=OMPID; Nc<GridN_Atom[Gc_AN]; Nc+=Nthrds){

	  GNc = Tmp_GridListAtom[Mc_AN][Nc];
          GRc = Tmp_CellListAtom[Mc_AN][Nc];

	  GN2N(GNc,N3);

	  n1 = N3[1] + Ngrid1*atv_ijk[GRc][1] - Min_Grid_Index[1];  
	  n2 = N3[2] + Ngrid2*atv_ijk[GRc][2] - Min_Grid_Index[2];  
	  n3 = N3[3] + Ngrid3*atv_ijk[GRc][3] - Min_Grid_Index[3];  

	  MGridListAtom[Mc_AN][Nc] = n1*Ng2*Ng3 + n2*Ng3 + n3;
	}

      } /* #pragma omp parallel */

    }

  }

  /* for PrintMemory */
  firsttime = 0;

  /****************************************************
                          Free
  ****************************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);

  for (Mc_AN=0; Mc_AN<=(Matomnum+MatomnumF); Mc_AN++){
    free(Tmp_CellListAtom[Mc_AN]);
    free(Tmp_GridListAtom[Mc_AN]);
  }
  free(Tmp_CellListAtom);
  free(Tmp_GridListAtom);

  if (alloc_first2==0 && estimate_switch!=0){

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(MGridListAtom[Mc_AN]);
    }
    free(MGridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(GridListAtom[Mc_AN]);
    }
    free(GridListAtom);

    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      free(CellListAtom[Mc_AN]);
    }
    free(CellListAtom);
  }

}




void Construct_MPI_Data_Structure_Grid()
{
  static int firsttime=1;
  int i,j,k,Mc_AN,Gc_AN,wan,n1,n2,n3;
  int min_n1,max_n1,min_n2,max_n2,N3[4];
  unsigned long long int AN,BN,CN,DN;
  unsigned long long int B_AB2,Bs,Be;
  unsigned long long int BN_AB,BN_CB,BN_CA,GN_B_AB;
  unsigned long long int GN,GNs,GR,n2D,N2D;
  unsigned long long int GN_AB,GN_CB,GN_CA;
  int size_Index_Snd_Grid_A2B;
  int size_Index_Rcv_Grid_A2B;
  int myid,numprocs,ID,IDS,IDR,tag=999;
  double Vec0,Vec1,coef,MinV,MaxV,rcut;
  double Cxyz[4],b[4],c[4],v[4];
  int alloc_first26;

  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /******************************************************
    find the smallest parallelepipedon which contains 
    atoms allocated to my ID under consideration of 
    cutoff radii of basis functions. 
  ******************************************************/
 
  for (k=1; k<=3; k++){

    if      (k==1){ i = 2; j = 3; }
    else if (k==2){ i = 3; j = 1; }
    else if (k==3){ i = 1; j = 2; }

    b[1] = tv[i][1];
    b[2] = tv[i][2];
    b[3] = tv[i][3];

    c[1] = tv[j][1];
    c[2] = tv[j][2];
    c[3] = tv[j][3];

    Cross_Product(b,c,v);
    coef = 1.0/sqrt(fabs( Dot_Product(v,v) ));
    
    v[1] = coef*v[1];
    v[2] = coef*v[2];
    v[3] = coef*v[3];

    MinV =  1.0e+10;
    MaxV = -1.0e+10;

    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
      Gc_AN = M2G[Mc_AN];
      wan  = WhatSpecies[Gc_AN];
      rcut = Spe_Atom_Cut1[wan];

      Cxyz[1] = Gxyz[Gc_AN][1] + rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] + rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] + rcut*v[3] - Grid_Origin[3];

      Vec0 = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      Cxyz[1] = Gxyz[Gc_AN][1] - rcut*v[1] - Grid_Origin[1];
      Cxyz[2] = Gxyz[Gc_AN][2] - rcut*v[2] - Grid_Origin[2];
      Cxyz[3] = Gxyz[Gc_AN][3] - rcut*v[3] - Grid_Origin[3];

      Vec1 = Dot_Product(Cxyz,rgtv[k])*0.5/PI;

      if (Vec0<MinV) MinV = Vec0;
      if (Vec1<MinV) MinV = Vec1;
      if (MaxV<Vec0) MaxV = Vec0;   
      if (MaxV<Vec1) MaxV = Vec1;
    }

    Min_Grid_Index[k] = (int)MinV;  /* buffer for GGA */
    Max_Grid_Index[k] = (int)MaxV;  /* buffer for GGA */

  } /* k */

  /******************************************************
    find the smallest parallelepipedon which contains 
    grids in the partition B_AB
    
    The parallelepipedon defines the partition D.
  ******************************************************/

  N2D = Ngrid1*Ngrid2;
  Bs = (myid*N2D+numprocs-1)/numprocs;
  Be = ((myid+1)*N2D+numprocs-1)/numprocs;
  
  min_n1 = 1000000;
  max_n1 =-1000000;
  min_n2 = 1000000;
  max_n2 =-1000000;

  for (B_AB2=Bs; B_AB2<Be; B_AB2++){

    n1 = B_AB2/Ngrid2;
    n2 = B_AB2 - n1*Ngrid2;

    if (n1<min_n1) min_n1 = n1;
    if (max_n1<n1) max_n1 = n1;
    if (n2<min_n2) min_n2 = n2;
    if (max_n2<n2) max_n2 = n2;
  }

  Min_Grid_Index_D[1] = min_n1 - 2;
  Max_Grid_Index_D[1] = max_n1 + 2;
  Min_Grid_Index_D[2] = min_n2 - 2;
  Max_Grid_Index_D[2] = max_n2 + 2;
  Min_Grid_Index_D[3] = -2;
  Max_Grid_Index_D[3] = (Ngrid3-1) + 2;

  /****************************************************************
      The partitions A to B

      construct the data structure for transfering rho_i from 
      the partitions A to B when rho is calculated 
      in the partition B using rho_i 
  ****************************************************************/

  /* find Num_Snd_Grid_A2B[ID] */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan  = WhatSpecies[Gc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];

      /* get process ID and increment Num_Snd_Grid_A2B */

      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = n2D*numprocs/N2D;
      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: Num_Snd_Grid_A2B */  

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
    MPI_Isend(&Num_Snd_Grid_A2B[IDS], 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
    MPI_Recv(&Num_Rcv_Grid_A2B[IDR], 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    MPI_Wait(&request,&stat);
  }

  /* allocation of arrays */  

  size_Index_Snd_Grid_A2B = 0;
  Index_Snd_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Snd_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Snd_Grid_A2B[ID]);
    size_Index_Snd_Grid_A2B += 3*Num_Snd_Grid_A2B[ID];
  }  
  
  size_Index_Rcv_Grid_A2B = 0; 
  Index_Rcv_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Index_Rcv_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Rcv_Grid_A2B[ID]);
    size_Index_Rcv_Grid_A2B += 3*Num_Rcv_Grid_A2B[ID]; 
  }  


  /* construct Index_Snd_Grid_A2B */  

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;

  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];
    wan  = WhatSpecies[Gc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GR = CellListAtom[Mc_AN][AN];

      /* get process ID and grid index (BN) for the partition B */

      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = n2D*numprocs/N2D;
      GN_B_AB = N3[1]*Ngrid2*Ngrid3 + N3[2]*Ngrid3 + N3[3];
      BN = GN_B_AB - ((ID*N2D+numprocs-1)/numprocs)*Ngrid3;

      /*
      printf("ABC1 myid=%2d ID=%2d\n",myid,ID);
      */

      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+0] = BN; 
      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+1] = Gc_AN; 
      Index_Snd_Grid_A2B[ID][3*Num_Snd_Grid_A2B[ID]+2] = GR; 

      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: Index_Snd_Grid_A2B */
      
  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    tag = 999;
   
    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Index_Snd_Grid_A2B[IDS][0], 3*Num_Snd_Grid_A2B[IDS], 
                 MPI_INT, IDS, tag, mpi_comm_level1, &request);
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Recv( &Index_Rcv_Grid_A2B[IDR][0], 3*Num_Rcv_Grid_A2B[IDR], 
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
    }

    if (Num_Snd_Grid_A2B[IDS]!=0)  MPI_Wait(&request,&stat);
  }

  /* count NN_A2B_S and NN_A2B_R */
  
  NN_A2B_S = 0;
  NN_A2B_R = 0;

  for (ID=0; ID<numprocs; ID++){
    if (Num_Snd_Grid_A2B[ID]!=0) NN_A2B_S++;
    if (Num_Rcv_Grid_A2B[ID]!=0) NN_A2B_R++;
  }

  N2D = Ngrid1*Ngrid2;
  My_NumGridB_AB =  (((myid+1)*N2D+numprocs-1)/numprocs)*Ngrid3
                  - ((myid*N2D+numprocs-1)/numprocs)*Ngrid3;

}


