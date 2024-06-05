/**********************************************************************
  Set_ddOrbitals_Grid_fuku.c:

   Set_ddOrbitals_Grid_fuku.c is a subroutine to calculate the value of basis
   functions on each grid point.

  Log of Set_ddOrbitals_Grid_fuku.c:

     23/Jan/2015  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "flpq.h"
#include "mpi.h"
#include <omp.h>




double Set_ddOrbitals_Grid()
{
  int i,j,k,n,Mc_AN,Gc_AN,Cwan,NO0,GNc,GRc;
  int Gh_AN,Mh_AN,Rnh,Hwan,NO1,Nog,h_AN;
  long int Nc;
  double time0;
  double x,y,z;
  double TStime,TEtime;
  double Cxyz[4];
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom,Etime_atom;
  int Norder_Chev=5000;       /* order of chebyshev expamsion to calculate its coefficients */

  MPI_Status stat;
  MPI_Request request;

  /* for OpenMP */
  int OMPID,Nthrds,Nprocs;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
//
//  Chev_Orbs = (double****)malloc(sizeof(double***)*SpeciesNum);
//  for (i=0; i<SpeciesNum; i++){
//    Chev_Orbs[i] = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      Chev_Orbs[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
//      for (k=0; k<List_YOUSO[24]; k++){
//        Chev_Orbs[i][j][k] = (double*)malloc(sizeof(double)*Norder_Chev);
//      }
//    }
//  }
//
//  Chev_dOrbs = (double****)malloc(sizeof(double***)*SpeciesNum);
//  for (i=0; i<SpeciesNum; i++){
//    Chev_dOrbs[i] = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      Chev_dOrbs[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
//      for (k=0; k<List_YOUSO[24]; k++){
//        Chev_dOrbs[i][j][k] = (double*)malloc(sizeof(double)*Norder_Chev);
//      }
//    }
//  }
//
//  Chev_ddOrbs = (double****)malloc(sizeof(double***)*SpeciesNum);
//  for (i=0; i<SpeciesNum; i++){
//    Chev_ddOrbs[i] = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      Chev_ddOrbs[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
//      for (k=0; k<List_YOUSO[24]; k++){
//        Chev_ddOrbs[i][j][k] = (double*)malloc(sizeof(double)*Norder_Chev);
//      }
//    }
//  }


  /*****************************************************
                Calculate orbitals on grids
  *****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];    
    Cwan = WhatSpecies[Gc_AN];

    if (Cnt_kind==0)  NO0 = Spe_Total_NO[Cwan];
    else              NO0 = Spe_Total_CNO[Cwan]; 

    /* calculate Chebyshev coefficients */
    //printf("start Get_Chebyshev_Coef\n");
    //Get_Chebyshev_Coef(Cwan, Norder_Chev, Chev_Orbs, Chev_dOrbs, Chev_ddOrbs);
    //printf("finish Get_Chebyshev_Coef\n");
    //printf("start Get_Chebyshev_Coef_XV\n");
    //Get_Chebyshev_Coef_XV(Cwan, Norder_Chev, Chev_Orbs, Chev_dOrbs, Chev_ddOrbs);
    //printf("finish Get_Chebyshev_Coef_XV\n");


//#pragma omp parallel shared(List_YOUSO,Orbs_Grid,dOrbs_Grid,ddOrbs_Grid,Cnt_kind,Gxyz,atv,CellListAtom,GridListAtom,GridN_Atom,Gc_AN,Cwan,Mc_AN,NO0,Chev_Orbs,Chev_dOrbs,Chev_ddOrbs) private(OMPID,Nthrds,Nprocs,Nc,GNc,GRc,Cxyz,x,y,z,i,j)
#pragma omp parallel shared(List_YOUSO,Orbs_Grid,dOrbs_Grid,ddOrbs_Grid,Cnt_kind,Gxyz,atv,CellListAtom,GridListAtom,GridN_Atom,Gc_AN,Cwan,Mc_AN,NO0) private(OMPID,Nthrds,Nprocs,Nc,GNc,GRc,Cxyz,x,y,z,i,j)
    {
      double **ddChi;
      double Cxyz0[4]; 
      int i,j;

      /* allocation of array */

      ddChi = (double**)malloc(sizeof(double*)*11);
      for (j=0; j<11; j++){
         ddChi[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      }

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();
      Nprocs = omp_get_num_procs();

      for (Nc=OMPID*GridN_Atom[Gc_AN]/Nthrds; Nc<(OMPID+1)*GridN_Atom[Gc_AN]/Nthrds; Nc++){

	GNc = GridListAtom[Mc_AN][Nc]; 
	GRc = CellListAtom[Mc_AN][Nc];

	Get_Grid_XYZ(GNc,Cxyz);
	x = Cxyz[1] + atv[GRc][1] - Gxyz[Gc_AN][1]; 
	y = Cxyz[2] + atv[GRc][2] - Gxyz[Gc_AN][2]; 
	z = Cxyz[3] + atv[GRc][3] - Gxyz[Gc_AN][3];

	if (Cnt_kind==0){
          Get_ddOrbitals(Cwan,x,y,z,ddChi);
          //Get_ddOrbitals_QuinticSpline(Cwan,x,y,z,ddChi);
          //Get_Chebyshev_Orbs(Cwan,x,y,z,ddChi);
          //Get_Chebyshev_Orbs_XV(Cwan,x,y,z,ddChi);
	}
	else{
          Get_Cnt_ddOrbitals(Mc_AN,x,y,z,ddChi);
	}

	for (i=0; i<NO0; i++){
	   Orbs_Grid[Mc_AN][Nc][i]    = (Type_Orbs_Grid)ddChi[0][i];/* AITUNE */
	  dOrbs_Grid[0][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[1][i];/* AITUNE */
	  dOrbs_Grid[1][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[2][i];/* AITUNE */
	  dOrbs_Grid[2][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[3][i];/* AITUNE */
	 ddOrbs_Grid[0][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[4][i];/* AITUNE */
  	 ddOrbs_Grid[1][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[5][i];/* AITUNE */
  	 ddOrbs_Grid[2][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[6][i];/* AITUNE */
  	 ddOrbs_Grid[3][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[7][i];/* AITUNE */
  	 ddOrbs_Grid[4][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[8][i];/* AITUNE */
  	 ddOrbs_Grid[5][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[9][i];/* AITUNE */
  	 ddOrbs_Grid[6][Mc_AN][Nc][i] = (Type_Orbs_Grid)ddChi[10][i];/* AITUNE */
	}

      } /* Nc */

      /* freeing of array */

      for (j=0; j<11; j++){
         free(ddChi[j]);
      }
      free(ddChi);

    } /* #pragma omp parallel */

  }

  /****************************************************
     Calculate Orbs_Grid_FNAN
  ****************************************************/

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];    

    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

      Gh_AN = natn[Gc_AN][h_AN];

      if (G2ID[Gh_AN]!=myid){

        Mh_AN = F_G2M[Gh_AN];
        Rnh = ncn[Gc_AN][h_AN];
        Hwan = WhatSpecies[Gh_AN];

        if (Cnt_kind==0)  NO1 = Spe_Total_NO[Hwan];
        else              NO1 = Spe_Total_CNO[Hwan];

        ///* calculate Chebyshev coefficients */
        //printf("start Get_Chebyshev_Coef in FNAN\n");
        //Get_Chebyshev_Coef(Hwan, Norder_Chev, Chev_Orbs);
        //printf("finish Get_Chebyshev_Coef in FNAN\n");

#pragma omp parallel shared(List_YOUSO,Orbs_Grid_FNAN,dOrbs_Grid_FNAN,ddOrbs_Grid_FNAN,NO1,Mh_AN,Hwan,Cnt_kind,Rnh,Gh_AN,Gxyz,atv,NumOLG,Mc_AN,h_AN,GListTAtoms1,GridListAtom,CellListAtom) private(OMPID,Nthrds,Nprocs,Nog,Nc,GNc,GRc,x,y,z,j)
        {

     double **ddChi;
	  double Cxyz0[4]; 

          /* allocation of arrays */

      ddChi = (double**)malloc(sizeof(double*)*11);
      for (j=0; j<11; j++){
         ddChi[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      }

	  /* get info. on OpenMP */ 

	  OMPID = omp_get_thread_num();
	  Nthrds = omp_get_num_threads();
	  Nprocs = omp_get_num_procs();

	  for (Nog=OMPID*NumOLG[Mc_AN][h_AN]/Nthrds; Nog<(OMPID+1)*NumOLG[Mc_AN][h_AN]/Nthrds; Nog++){

	    Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	    GNc = GridListAtom[Mc_AN][Nc];
	    GRc = CellListAtom[Mc_AN][Nc]; 

	    Get_Grid_XYZ(GNc,Cxyz0);

	    x = Cxyz0[1] + atv[GRc][1] - Gxyz[Gh_AN][1] - atv[Rnh][1];
	    y = Cxyz0[2] + atv[GRc][2] - Gxyz[Gh_AN][2] - atv[Rnh][2];
	    z = Cxyz0[3] + atv[GRc][3] - Gxyz[Gh_AN][3] - atv[Rnh][3];

	    if (Cnt_kind==0){
              Get_ddOrbitals(Hwan,x,y,z,ddChi);
              //Get_ddOrbitals_QuinticSpline(Hwan,x,y,z,ddChi);
              //Get_Chebyshev_Orbs(Hwan,x,y,z,ddChi);
	    } 
	    else{
              Get_Cnt_ddOrbitals(Mh_AN,x,y,z,ddChi);
	    }

	    for (j=0; j<NO1; j++){
	      Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j]     = (Type_Orbs_Grid)ddChi[0][j];/* AITUNE */
	      dOrbs_Grid_FNAN[0][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[1][j];/* AITUNE */
	      dOrbs_Grid_FNAN[1][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[2][j];/* AITUNE */
	      dOrbs_Grid_FNAN[2][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[3][j];/* AITUNE */
	     ddOrbs_Grid_FNAN[0][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[4][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[1][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[5][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[2][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[6][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[3][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[7][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[4][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[8][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[5][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[9][j];/* AITUNE */
  	     ddOrbs_Grid_FNAN[6][Mc_AN][h_AN][Nog][j] = (Type_Orbs_Grid)ddChi[10][j];/* AITUNE */
	    }

	  } /* Nog */

          /* freeing of arrays */
          for (j=0; j<11; j++){
             free(ddChi[j]);
          }
          free(ddChi);
        } 
      }
    }
  }


//  for (i=0; i<SpeciesNum; i++){
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      for (k=0; k<List_YOUSO[24]; k++){
//        free(Chev_Orbs[i][j][k]);
//      }
//      free(Chev_Orbs[i][j]);
//    }
//    free(Chev_Orbs[i]);
//  }
//  free(Chev_Orbs);
//
//  for (i=0; i<SpeciesNum; i++){
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      for (k=0; k<List_YOUSO[24]; k++){
//        free(Chev_dOrbs[i][j][k]);
//      }
//      free(Chev_dOrbs[i][j]);
//    }
//    free(Chev_dOrbs[i]);
//  }
//  free(Chev_dOrbs);
//
//  for (i=0; i<SpeciesNum; i++){
//    for (j=0; j<(List_YOUSO[25]+1); j++){
//      for (k=0; k<List_YOUSO[24]; k++){
//        free(Chev_ddOrbs[i][j][k]);
//      }
//      free(Chev_ddOrbs[i][j]);
//    }
//    free(Chev_ddOrbs[i]);
//  }
//  free(Chev_ddOrbs);

  return time0;
}
