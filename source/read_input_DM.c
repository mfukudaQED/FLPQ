/**********************************************************************
  read_input_DM.c:

     read_input_DM.c is a subroutine to read DM from
     OpenMX_LPQ_DM_*.bin. 
     

  Log of read_input_DM.c:

     21/May/2019  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"

void read_input_DM( int kloop, int mu )
{

  int i,j,k;
  int numprocs,myid,tag=999,ID;
  FILE *fp;
  char fname[500];
  double *****DM_full;
  double *****iDM_full;


  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);


  {
    int spin,Mc_AN,Gc_AN,wan1,TNO1;
    int h_AN,Gh_AN,wan2,TNO2;

    DM_full = (double*****)malloc(sizeof(double****)*(SpinP_switch+1));
    for (spin=0; spin<=SpinP_switch; spin++){

      DM_full[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
      for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
        if(Gc_AN==0){
          TNO1=1;
          FNAN[0]=0;
        }
        else{
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
        }
        DM_full[spin][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          DM_full[spin][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

          if (Gc_AN==0){ 
            TNO2 = 1;
          }
          else{ 
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
          }
          for (i=0; i<TNO1; i++){
            DM_full[spin][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
          }
        }
      }
    }

    for (spin=0; spin<=SpinP_switch; spin++){
      /*** open input DM ***/
      //open_input_DM(fp,spin,kloop,mu);
      if(flag_DMmu_range==1){
        sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DMmu_%d.bin",spin);
        if ((fp = fopen(fname,"rb")) == NULL){
          printf("cannot open %s \n", fname);
          exit(0);
        }
      }

      else if((flag_DMmu==0)&&(flag_DMk==0)){
        sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d.bin",spin);
        if ((fp = fopen(fname,"rb")) == NULL){
          printf("cannot open %s \n", fname);
          exit(0);
        }
      }

      else if((flag_DMmu==1)&&(flag_DMk==0)){
        if(spin==Num_DM_spin){
          sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d_%08d.bin",spin,mu);
          if ((fp = fopen(fname,"rb")) == NULL){
            printf("cannot open %s \n", fname);
            exit(0);
          }
        }
      }

      else if((flag_DMmu==1)&&(flag_DMk==1)){
        if(spin==Num_DM_spin){
          sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d_%08d_%08d.bin",spin,kloop,mu);
          //printf("read %s\n",fname);
          if ((fp = fopen(fname,"rb")) == NULL){
            printf("cannot open %s \n", fname);
            exit(0);
          }
        }
      }

      else {
        printf("error flag_DMmu, flag_DMk\n");
        exit(0);
      }

      //printf("read %s\n",fname);
      /*** ***/

      if((flag_DMmu==1)&&(spin!=Num_DM_spin)){

        for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
          for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
            for (i=0; i<TNO1; i++){
              for (j=0; j<TNO2; j++){
                DM_full[spin][Gc_AN][h_AN][i][j] = 0.0;
              }
            }
          }
        }
      }

      else{
    
        for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
          for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
            for (i=0; i<TNO1; i++){
              fread(DM_full[spin][Gc_AN][h_AN][i],sizeof(double),TNO2,fp);
              //for (j=0; j<TNO2; j++){
              //  printf("%lf\n",DM_full[spin][Gc_AN][h_AN][i][j]);
              //}
            }
          }
        }
        fclose(fp);
      }

    }


  }

  /* allocate Density Matrix */
  {
    int spin,Mc_AN,Gc_AN,Cwan,tno0;
    int h_AN,Gh_AN,Hwan,tno1;

    DM = (double******)malloc(sizeof(double*****)*1); 
    DM[0] = (double*****)malloc(sizeof(double****)*(SpinP_switch+1)); 
    for (spin=0; spin<=SpinP_switch; spin++){
      DM[0][spin] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
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

	DM[0][spin][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;  
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_NO[Hwan];
	  } 

	  DM[0][spin][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
	  for (i=0; i<tno0; i++){
	    DM[0][spin][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
  	    for (j=0; j<tno1; j++) DM[0][spin][Mc_AN][h_AN][i][j] = DM_full[spin][Gc_AN][h_AN][i][j]; 
	  }
	}
      }
    }
  }

  {
    int spin,Mc_AN,Gc_AN,wan1,TNO1;
    int h_AN,Gh_AN,wan2,TNO2;

    for (spin=0; spin<=SpinP_switch; spin++){

      for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
        if(Gc_AN==0){
          TNO1=1;
          FNAN[0]=0;
        }
        else{
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
        }
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          for (i=0; i<TNO1; i++){
            free(DM_full[spin][Gc_AN][h_AN][i]);
          }
          free(DM_full[spin][Gc_AN][h_AN]);
        }
        free(DM_full[spin][Gc_AN]);
      }
      free(DM_full[spin]);
    }
    free(DM_full);
  }

  ///* read Density Matrix */
  //{
  //  int spin,MA_AN,GA_AN,wanA,tnoA;
  //  int LB_AN,GB_AN,wanB,tnoB;

  //  for (spin=0; spin<=SpinP_switch; spin++) {
  //    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++) {
  //      GA_AN = M2G[MA_AN];
  //      wanA = WhatSpecies[GA_AN];
  //      tnoA = Spe_Total_CNO[wanA];
  //      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
  //        GB_AN = natn[GA_AN][LB_AN];
  //        wanB = WhatSpecies[GB_AN];
  //        tnoB = Spe_Total_CNO[wanB];
  //        for (i=0; i<tnoA; i++){
  //          fread(DM[0][spin][MA_AN][LB_AN][i],sizeof(double),tnoB,fp);
  //        }
  //      }
  //    }
  //  }
  //}



  if ((Solver==4)||(SpinP_switch>1)){

  {
    int spin,Mc_AN,Gc_AN,wan1,TNO1;
    int h_AN,Gh_AN,wan2,TNO2;

    iDM_full = (double*****)malloc(sizeof(double****)*(2));
    for (spin=0; spin<2; spin++){

      iDM_full[spin] = (double****)malloc(sizeof(double***)*(atomnum+1));
      for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
        if(Gc_AN==0){
          TNO1=1;
          FNAN[0]=0;
        }
        else{
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
        }
        iDM_full[spin][Gc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          iDM_full[spin][Gc_AN][h_AN] = (double**)malloc(sizeof(double*)*TNO1);

          if (Gc_AN==0){ 
            TNO2 = 1;
          }
          else{ 
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
          }
          for (i=0; i<TNO1; i++){
            iDM_full[spin][Gc_AN][h_AN][i] = (double*)malloc(sizeof(double)*TNO2);
          }
        }
      }
    }

    for (spin=0; spin<2; spin++){
      /*** open input iDM ***/
      if(flag_DMmu_range==1){
        sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DMmu_%d.bin",spin+SpinP_switch+1);
        if ((fp = fopen(fname,"rb")) == NULL){
          printf("cannot open %s \n", fname);
          exit(0);
        }
      }

      else if((flag_DMmu==0)&&(flag_DMk==0)){
        sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d.bin",spin+SpinP_switch+1);
        if ((fp = fopen(fname,"rb")) == NULL){
          printf("cannot open %s \n", fname);
          exit(0);
        }
      }

      else if((flag_DMmu==1)&&(flag_DMk==0)){
        if(spin==Num_DM_spin){
          sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d_%08d.bin",spin+SpinP_switch+1,mu);
          if ((fp = fopen(fname,"rb")) == NULL){
            printf("cannot open %s \n", fname);
            exit(0);
          }
        }
      }

      else if((flag_DMmu==1)&&(flag_DMk==1)){
        if(spin==Num_DM_spin){
          sprintf(fname,"OpenMX_LPQ_DM/OpenMX_LPQ_DM_%d_%08d_%08d.bin",spin+SpinP_switch+1,kloop,mu);
          if ((fp = fopen(fname,"rb")) == NULL){
            printf("cannot open %s \n", fname);
            exit(0);
          }
        }
      }

      else {
        printf("error flag_DMmu, flag_DMk\n");
        exit(0);
      }
      printf("read %s\n",fname);
      /*** ***/

      if((flag_DMmu==1)&&(spin!=Num_DM_spin)){

        for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
          for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
            for (i=0; i<TNO1; i++){
              fread(iDM_full[spin][Gc_AN][h_AN][i],sizeof(double),TNO2,fp);
              for (j=0; j<TNO2; j++){
                iDM_full[spin][Gc_AN][h_AN][i][j] = 0.0;
              }
            }
          }
        }
      }
      else{
      
        for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
          for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
            Gh_AN = natn[Gc_AN][h_AN];
            wan2 = WhatSpecies[Gh_AN];
            TNO2 = Spe_Total_CNO[wan2];
            for (i=0; i<TNO1; i++){
              fread(iDM_full[spin][Gc_AN][h_AN][i],sizeof(double),TNO2,fp);
              //for (j=0; j<TNO2; j++){
              //  printf("%18.8E\n",iDM_full[spin][Gc_AN][h_AN][i][j]);
              //}
            }
          }
        }
        fclose(fp);
      }
      
    }

  }

  {
    int spin,Mc_AN,Gc_AN,Cwan,tno0;
    int h_AN,Gh_AN,Hwan,tno1;

    iDM = (double******)malloc(sizeof(double*****)*1);
    iDM[0] = (double*****)malloc(sizeof(double****)*2); 
    for (spin=0; spin<2; spin++){
      iDM[0][spin] = (double****)malloc(sizeof(double***)*(Matomnum+1)); 
      FNAN[0] = 0;
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

        iDM[0][spin][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1)); 
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

          if (Mc_AN==0){
            tno1 = 1;  
          }
          else{
            Gh_AN = natn[Gc_AN][h_AN];
            Hwan = WhatSpecies[Gh_AN];
            tno1 = Spe_Total_NO[Hwan];
          } 

          iDM[0][spin][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0); 
          for (i=0; i<tno0; i++){
            iDM[0][spin][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1); 
            for (j=0; j<tno1; j++)  iDM[0][spin][Mc_AN][h_AN][i][j] = iDM_full[spin][Gc_AN][h_AN][i][j]; 
          }
        }
      }
    }
  }

  {
    int spin,Mc_AN,Gc_AN,wan1,TNO1;
    int h_AN,Gh_AN,wan2,TNO2;

    for (spin=0; spin<2; spin++){

      for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
        if(Gc_AN==0){
          TNO1=1;
          FNAN[0]=0;
        }
        else{
          wan1 = WhatSpecies[Gc_AN];
          TNO1 = Spe_Total_CNO[wan1];
        }
        for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
          for (i=0; i<TNO1; i++){
            free(iDM_full[spin][Gc_AN][h_AN][i]);
          }
          free(iDM_full[spin][Gc_AN][h_AN]);
        }
        free(iDM_full[spin][Gc_AN]);
      }
      free(iDM_full[spin]);
    }
    free(iDM_full);
  }

  ///* read imaginary parts of Density Matrix */
  //{
  //  int spin,MA_AN,GA_AN,wanA,tnoA;
  //  int LB_AN,GB_AN,wanB,tnoB;

  //  for (spin=0; spin<2; spin++) {
  //    for (MA_AN=1; MA_AN<=Matomnum; MA_AN++) {
  //      GA_AN = M2G[MA_AN];
  //      wanA = WhatSpecies[GA_AN];
  //      tnoA = Spe_Total_CNO[wanA];
  //      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
  //        GB_AN = natn[GA_AN][LB_AN];
  //        wanB = WhatSpecies[GB_AN];
  //        tnoB = Spe_Total_CNO[wanB];
  //        for (i=0; i<tnoA; i++){
  //          fread(iDM[0][spin][MA_AN][LB_AN][i],sizeof(double),tnoB,fp);
  //        }
  //      }
  //    }
  //  }
  //}

  }/* if ((Solver==4)||(SpinP_switch>1)) */

}
