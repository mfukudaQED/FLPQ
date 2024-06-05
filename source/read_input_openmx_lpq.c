/**********************************************************************
  read_input_openmx_lpq.c:

     read_input_openmx_lpq.c is a subroutine to read an input file from
     OpenMX_LPQ.bin. 
     Set_Allocate_Atom2CPU and Set_Inf_SndRcv are executed in read_input_openmx_lpq.
     

  Log of read_input_openmx_lpq.c:

     28/Jun/2016  Released by M.Fukuda
     20/May/2019  Modified by M.Fukuda for MPI

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"

void read_input_openmx_lpq()
{

  int i,j,k;
  int numprocs,myid,tag=999,ID;
  char fname[500];
  FILE *fp;
  double ***CntCoes_full;

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

/**********************************/
/* input data from OpenMX_LPQ.dat */
/**********************************/
  sprintf(fname,"OpenMX_LPQ.bin");

  if ((fp = fopen(fname,"rb")) == NULL){
    printf("cannot open %s \n", fname);
    exit(0);
  }

  fread(&SpinP_switch      ,sizeof(int),1,fp);
  fread(&Solver            ,sizeof(int),1,fp);
  fread(&Cnt_switch        ,sizeof(int),1,fp);
  fread(&Cnt_kind          ,sizeof(int),1,fp);
  fread(&Calc_CntOrbital_ON,sizeof(int),1,fp);
  fread(&MD_switch         ,sizeof(int),1,fp);
  fread(&level_stdout      ,sizeof(int),1,fp);
  fread(&ESM_switch        ,sizeof(int),1,fp);
  fread(&TCpyCell          ,sizeof(int),1,fp);
  fread(&CpyCell           ,sizeof(int),1,fp);
  fread(&Ngrid1,sizeof(int),1,fp);
  fread(&Ngrid2,sizeof(int),1,fp);
  fread(&Ngrid3,sizeof(int),1,fp);
  fread(List_YOUSO,sizeof(int),NYOUSO,fp);
  fread(&SpeciesNum,sizeof(int),1,fp);
  fread(&Max_FSNAN ,sizeof(int),1,fp);
  //fread(&Max_NumOLG,sizeof(int),1,fp);
  //fread(&Matomnum      ,sizeof(int),1,fp);
  //fread(&MatomnumF     ,sizeof(int),1,fp);
  fread(&atomnum       ,sizeof(int),1,fp);
  //fread(&NN_A2B_S      ,sizeof(int),1,fp);
  //fread(&NN_A2B_R      ,sizeof(int),1,fp);
  fread(&ScaleSize,sizeof(double),1,fp);
  fread(tv,sizeof(double),16,fp);
  fread(rtv,sizeof(double),16,fp);
  //for (i=0;i<4;i++){
  //  for (j=0;j<4;j++){
  //    fread(tv[i][j],sizeof(double),1,fp);
  //  }
  //}
  //for (i=0;i<4;i++){
  //  for (j=0;j<4;j++){
  //    fread(rtv[i][j],sizeof(double),1,fp);
  //  }
  //}

  Gxyz = (double**)malloc(sizeof(double*)*(atomnum+1));
  for (i=0; i<(atomnum+1); i++){
    Gxyz[i] = (double*)malloc(sizeof(double)*YOUSO26);
    fread(Gxyz[i],sizeof(double),YOUSO26,fp);
  }

  G2ID = (int*)malloc(sizeof(int)*(atomnum+1));
  //fread(G2ID,sizeof(int),atomnum+1,fp);
  
  /* Set Matomnum, M2G, G2ID */
  //printf("start Set_Allocate_Atom2CPU\n");
  Set_Allocate_Atom2CPU(0);
  //printf("end Set_Allocate_Atom2CPU\n");

  Spe_Total_NO = (int*)malloc(sizeof(int)*SpeciesNum);
  Spe_Total_CNO = (int*)malloc(sizeof(int)*SpeciesNum);
  Spe_Num_Mesh_PAO = (int*)malloc(sizeof(int)*SpeciesNum);
  Spe_MaxL_Basis = (int*)malloc(sizeof(int)*SpeciesNum);

  fread(Spe_Total_NO,    sizeof(int),SpeciesNum,fp);
  fread(Spe_Total_CNO,   sizeof(int),SpeciesNum,fp);
  fread(Spe_Num_Mesh_PAO,sizeof(int),SpeciesNum,fp);
  fread(Spe_MaxL_Basis,  sizeof(int),SpeciesNum,fp);


  Spe_Num_Basis = (int**)malloc(sizeof(int*)*SpeciesNum);
  for (i=0; i<SpeciesNum; i++){
    Spe_Num_Basis[i] = (int*)malloc(sizeof(int)*(Supported_MaxL+1));
    fread(Spe_Num_Basis[i],sizeof(int),Supported_MaxL+1,fp);
  }  
   
  Spe_Num_CBasis = (int**)malloc(sizeof(int*)*SpeciesNum);
  for (i=0; i<SpeciesNum; i++){
    Spe_Num_CBasis[i] = (int*)malloc(sizeof(int)*(Supported_MaxL+1));
    fread(Spe_Num_CBasis[i],sizeof(int),Supported_MaxL+1,fp);
  }  
  
  Spe_PAO_RV = (double**)malloc(sizeof(double*)*List_YOUSO[18]);
  for (i=0; i<List_YOUSO[18]; i++){
    Spe_PAO_RV[i] = (double*)malloc(sizeof(double)*List_YOUSO[21]);
    fread(Spe_PAO_RV[i],sizeof(double),List_YOUSO[21],fp);
  }

  Spe_PAO_RWF = (double****)malloc(sizeof(double***)*List_YOUSO[18]);
  for (i=0; i<List_YOUSO[18]; i++){
    Spe_PAO_RWF[i] = (double***)malloc(sizeof(double**)*(List_YOUSO[25]+1));
    for (j=0; j<=List_YOUSO[25]; j++){
      Spe_PAO_RWF[i][j] = (double**)malloc(sizeof(double*)*List_YOUSO[24]);
      for (k=0; k<List_YOUSO[24]; k++){
        Spe_PAO_RWF[i][j][k] = (double*)malloc(sizeof(double)*List_YOUSO[21]);
        fread(Spe_PAO_RWF[i][j][k],sizeof(double),List_YOUSO[21],fp);
      }
    }
  }

  Spe_Atom_Cut1 = (double*)malloc(sizeof(double)*SpeciesNum);
  fread(Spe_Atom_Cut1, sizeof(double),SpeciesNum,fp);
  //for (i=0; i<SpeciesNum; i++){
  //  printf("Spe_Atom_Cut1 %d, %lf\n",i,Spe_Atom_Cut1[i]);
  //}

  Spe_Specified_Num = (int**)malloc(sizeof(int*)*List_YOUSO[18]); 
  for (i=0; i<SpeciesNum; i++){
    Spe_Specified_Num[i] = (int*)malloc(sizeof(int)*Spe_Total_NO[i]); 
    fread(Spe_Specified_Num[i],sizeof(int),Spe_Total_NO[i],fp);
  }

  Spe_Trans_Orbital = (int***)malloc(sizeof(int**)*List_YOUSO[18]); 
  for (i=0; i<SpeciesNum; i++){
    Spe_Trans_Orbital[i] = (int**)malloc(sizeof(int*)*Spe_Total_NO[i]); 
    for (j=0; j<Spe_Total_NO[i]; j++){
      Spe_Trans_Orbital[i][j] = (int*)malloc(sizeof(int)*List_YOUSO[24]); 
      fread(Spe_Trans_Orbital[i][j],sizeof(int),List_YOUSO[24],fp);
    }
  }

  /* for Print_CubeTitle */
  Spe_WhatAtom = (int*)malloc(sizeof(int)*SpeciesNum);
  fread(Spe_WhatAtom, sizeof(int),SpeciesNum,fp);

  /* for Print_CubeTitle */
  Spe_Core_Charge = (double*)malloc(sizeof(double)*SpeciesNum);
  fread(Spe_Core_Charge, sizeof(double),SpeciesNum,fp);

  /* for Print_CubeTitle */
  InitN_USpin = (double*)malloc(sizeof(double)*(atomnum+1));
  InitN_DSpin = (double*)malloc(sizeof(double)*(atomnum+1));
  fread(InitN_USpin, sizeof(double),(atomnum+1),fp);
  fread(InitN_DSpin, sizeof(double),(atomnum+1),fp);
  
  //M2G = (int*)malloc(sizeof(int)*(Matomnum+2));
  //fread(M2G,sizeof(int),Matomnum+2,fp);

  //F_M2G = (int*)malloc(sizeof(int)*(Matomnum+MatomnumF+1));
  //fread(F_M2G,sizeof(int),Matomnum+MatomnumF+1,fp);
  
  //F_G2M = (int*)malloc(sizeof(int)*(atomnum+1));
  //fread(F_G2M,sizeof(int),atomnum+1,fp);

  WhatSpecies = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(WhatSpecies,sizeof(int),atomnum+1,fp);
  
  GridN_Atom = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(GridN_Atom,sizeof(int),atomnum+1,fp);
  
  FNAN = (int*)malloc(sizeof(int)*(atomnum+1));
  fread(FNAN,sizeof(int),atomnum+1,fp);

  ncn  = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (i=0; i<=atomnum; i++){
    ncn[i] = (int*)malloc(sizeof(int)*((int)(Max_FSNAN*ScaleSize)+1));
    fread(ncn[i],sizeof(int),(int)(Max_FSNAN*ScaleSize)+1,fp);
  }

  F_Snd_Num = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(F_Snd_Num,sizeof(int),Num_Procs,fp);

  //S_Snd_Num = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(S_Snd_Num,sizeof(int),Num_Procs,fp);

  F_Rcv_Num = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(F_Rcv_Num,sizeof(int),Num_Procs,fp);

  //S_Rcv_Num = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(S_Rcv_Num,sizeof(int),Num_Procs,fp);

  F_TopMAN = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(F_TopMAN,sizeof(int),Num_Procs,fp);

  //Snd_MAN = (int**)malloc(sizeof(int*)*numprocs);
  //for (i=0; i<numprocs; i++){
  //  Snd_MAN[i] = (int*)malloc(sizeof(int)*(F_Snd_Num[i]+1));
  //  //fread(Snd_MAN[i],sizeof(int),(F_Snd_Num[i]+S_Snd_Num[i]+1),fp);
  //}

  //Snd_GAN = (int**)malloc(sizeof(int*)*numprocs);
  //for (i=0; i<numprocs; i++){
  //  Snd_GAN[i] = (int*)malloc(sizeof(int)*(F_Snd_Num[i]+1));
  //  //fread(Snd_GAN[i],sizeof(int),(F_Snd_Num[i]+S_Snd_Num[i]+1),fp);
  //}

  //Rcv_GAN = (int**)malloc(sizeof(int*)*numprocs);
  //for (i=0; i<numprocs; i++){
  //  Rcv_GAN[i] = (int*)malloc(sizeof(int)*(F_Rcv_Num[i]+S_Rcv_Num[i]+1));
  //  //fread(Rcv_GAN[i],sizeof(int),(F_Rcv_Num[i]+S_Rcv_Num[i]+1),fp);
  //}

  Num_Snd_Grid_A2B = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(Num_Snd_Grid_A2B,sizeof(int),Num_Procs,fp);

  Num_Rcv_Grid_A2B = (int*)malloc(sizeof(int)*Num_Procs);
  //fread(Num_Rcv_Grid_A2B,sizeof(int),Num_Procs,fp);

  //Index_Snd_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  //for (ID=0; ID<numprocs; ID++){
  //  Index_Snd_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Snd_Grid_A2B[ID]);
  //  fread(Index_Snd_Grid_A2B[ID],sizeof(int),3*Num_Snd_Grid_A2B[ID],fp);
  //}  
  //
  //Index_Rcv_Grid_A2B = (int**)malloc(sizeof(int*)*numprocs);
  //for (ID=0; ID<numprocs; ID++){
  //  Index_Rcv_Grid_A2B[ID] = (int*)malloc(sizeof(int)*3*Num_Rcv_Grid_A2B[ID]);
  //  fread(Index_Rcv_Grid_A2B[ID],sizeof(int),3*Num_Rcv_Grid_A2B[ID],fp);
  //}  

  {
  int n,TN;
  TN = (2*CpyCell+1)*(2*CpyCell+1)*(2*CpyCell+1) - 1;

    atv = (double**)malloc(sizeof(double*)*(TN+1));
    for (i=0; i<(TN+1); i++){
      atv[i] = (double*)malloc(sizeof(double)*4);
      fread(atv[i],sizeof(double),4,fp);
    }

    n = 2*CpyCell + 4;
    ratv = (int***)malloc(sizeof(int**)*n);
    for (i=0; i<n; i++){
      ratv[i] = (int**)malloc(sizeof(int*)*n);
      for (j=0; j<n; j++){
	ratv[i][j] = (int*)malloc(sizeof(int)*n);
        fread(ratv[i][j],sizeof(int),n,fp);
      }
    }

    atv_ijk = (int**)malloc(sizeof(int*)*(TN+1));
    for (i=0; i<(TN+1); i++){
      atv_ijk[i] = (int*)malloc(sizeof(int)*4);
      fread(atv_ijk[i],sizeof(int),4,fp);
    }
  }

  natn = (int**)malloc(sizeof(int*)*(atomnum+1));
  for (i=0; i<=atomnum; i++){
    natn[i] = (int*)malloc(sizeof(int)*((int)(Max_FSNAN*ScaleSize)+1));
    fread(natn[i],sizeof(int),((int)(Max_FSNAN*ScaleSize)+1),fp);
  }

  fread(&flag_energy_range_DM,sizeof(int),1,fp);
  if(flag_energy_range_DM==1){
    fread(DM_energy_range,sizeof(double),2,fp);
  }


   /*  NumOLG is setted in Set_Grid.c */
  {
    int Mc_AN,Gc_AN;
    NumOLG = (int**)malloc(sizeof(int*)*(Matomnum+1));
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
      if (Mc_AN==0) Gc_AN = 0;
      else          Gc_AN = M2G[Mc_AN];
      NumOLG[Mc_AN] = (int*)malloc(sizeof(int)*(FNAN[Gc_AN]+1));
    }
  }

  /*  Cell_Gxyz is setted in Set_Grid.c */
  //Cell_Gxyz = (double**)malloc(sizeof(double*)*(atomnum+1));
  //for (i=0; i<(atomnum+1); i++){
  //  Cell_Gxyz[i] = (double*)malloc(sizeof(double)*4);
  //}      
 

  /*****************************************************
                 Set_Inf_SndRcv
  *****************************************************/
         Set_Inf_SndRcv();


  /* CntCoes */
  if (Cnt_switch==1){
    int num;
    int L0,Mul0,M0,al,p;
    int spin,Mc_AN,Gc_AN,wan;

    CntCoes_full = (double***)malloc(sizeof(double**)*(atomnum+1));
    for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
      CntCoes_full[Gc_AN] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
        CntCoes_full[Gc_AN][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
      }
    }
      
    for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){

      wan = WhatSpecies[Gc_AN];

      al = -1;
      for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
        for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
          for (M0=0; M0<=2*L0; M0++){
            al++;
            fread(CntCoes_full[Gc_AN][al],sizeof(double),Spe_Specified_Num[wan][al],fp);
          }
        }
      }
    }

    CntCoes = (double***)malloc(sizeof(double**)*(Matomnum+MatomnumF+1));
    for (i=0; i<=(Matomnum+MatomnumF); i++){
      CntCoes[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]);
      for (j=0; j<List_YOUSO[7]; j++){
        CntCoes[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
      }

      Gc_AN = F_M2G[i];
      wan = WhatSpecies[Gc_AN];

      al = -1;
      num = 0;
      for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
        for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
          for (M0=0; M0<=2*L0; M0++){
            al++;
            for (p=0; p<Spe_Specified_Num[wan][al]; p++){
              CntCoes[i][al][p] = CntCoes_full[Gc_AN][al][p];
            }
          }
        }
      }
    }

    for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
      for (j=0; j<List_YOUSO[7]; j++){
        free(CntCoes_full[Gc_AN][j]);
      }
      free(CntCoes_full[Gc_AN]);
    }
    free(CntCoes_full);
     
  }


  fclose(fp);

}

