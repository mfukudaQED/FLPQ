/**********************************************************************
  Set_dd2Density_Grid_fuku.c:

     Set_dd2Density_Grid_fuku.c is a subroutine to calculate local physical quantities 
     on grid by one-particle wave functions.

  Log of Set_dd2Density_Grid_fuku.c:

     28/Dec/2014  Released by M.Fukuda
     01/Jul/2016  Modified by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "flpq.h"
#include "mpi.h"
#include <omp.h>

#define  measure_time   0



double Set_dd2Density_Grid(double *****CDM)
{
  static int firsttime=1;
  int al,L0,Mul0,M0,p,size1,size2;
  int Gc_AN,Mc_AN,Mh_AN,LN,AN,BN,CN;
  int n1,n2,n3,k1,k2,k3,N3[4];
  int Cwan,NO0,NO1,Rn,N,Hwan,i,j,k,n;
  int NN_S,NN_R;
  unsigned long long int N2D,n2D,GN; 
  int Max_Size,My_Max;
  int size_Tmp_dDen_Grid;
  int size_Den_Snd_Grid_A2B;
  int size_Den_Rcv_Grid_A2B;
  int h_AN,Gh_AN,Rnh,spin,Nc,GRc,Nh,Nog;
  int Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3;

  double threshold;
  double tmp0,tmp1,sk1,sk2,sk3,tot_den,sum,sum2;
  double tmpx,tmpy,tmpz;
  double tmp0_0,tmp0_1,tmp0_2,tmp0_3;
  double sum_0,sum_1,sum_2,sum_3;
  double sum2_0,sum2_1,sum2_2,sum2_3;
  double d1,d2,d3,cop,sip,sit,cot;
  double x,y,z,Cxyz[4];
  double TStime,TEtime;
  double ****Tmp_dDen_Grid;
  double **Den_Snd_Grid_A2B;
  double **Den_Rcv_Grid_A2B;
  double *tmp_array,*tmpx_array,*tmpy_array,*tmpz_array;
  double *tmp_array2,*tmpx_array2,*tmpy_array2,*tmpz_array2;
  double *orbs0,*orbs1,**dorbs0,**dorbs1;
  double *orbs0_0,*orbs0_1,*orbs0_2,*orbs0_3;
  double *orbs1_0,*orbs1_1,*orbs1_2,*orbs1_3;
  double **dorbs0_0,**dorbs0_1,**dorbs0_2,**dorbs0_3;
  double **dorbs1_0,**dorbs1_1,**dorbs1_2,**dorbs1_3;
  double ***tmp_CDM;
  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;
  double time0,time1,time2;
  int SpinP_Grid_size;

  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  /* for OpenMP */
  int OMPID,Nthrds;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (SpinP_switch==3){
       SpinP_Grid_size = SpinP_switch + 2;
    }
  else{
       SpinP_Grid_size = SpinP_switch;
    }
  
//  if (SpinP_switch==3){
//      printf("Start Set_dDensity_Grid.\n");
//    }
//  else{
//      printf("Set_dDensity_Grid is unavailable for SpinP_switch<=2.\n");
//      exit(1);
//    }

  /* allocation of arrays */

  size_Tmp_dDen_Grid = 0;
  Tmp_dDen_Grid = (double****)malloc(sizeof(double***)*7); 
  for (k=0; k<7; k++){
     Tmp_dDen_Grid[k] = (double***)malloc(sizeof(double**)*(SpinP_Grid_size+1)); 
     for (i=0; i<(SpinP_Grid_size+1); i++){
       Tmp_dDen_Grid[k][i] = (double**)malloc(sizeof(double*)*(Matomnum+1)); 
       Tmp_dDen_Grid[k][i][0] = (double*)malloc(sizeof(double)*1); 
       for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
         Gc_AN = F_M2G[Mc_AN];
         Tmp_dDen_Grid[k][i][Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
        
         /* AITUNE */
         for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
            Tmp_dDen_Grid[k][i][Mc_AN][Nc] = 0.0;
         }
         
         size_Tmp_dDen_Grid += GridN_Atom[Gc_AN];
       }
     }
   }

  size_Den_Snd_Grid_A2B = 0; 
  Den_Snd_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Snd_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Snd_Grid_A2B[ID]*(SpinP_Grid_size+1)*7);
    size_Den_Snd_Grid_A2B += Num_Snd_Grid_A2B[ID]*(SpinP_Grid_size+1)*7;
  }  

  size_Den_Rcv_Grid_A2B = 0;   
  Den_Rcv_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Rcv_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Rcv_Grid_A2B[ID]*(SpinP_Grid_size+1)*7);
    size_Den_Rcv_Grid_A2B += Num_Rcv_Grid_A2B[ID]*(SpinP_Grid_size+1)*7;   
  }

  /****************************************************
                when orbital optimization
  ****************************************************/

  if (Cnt_kind==0 && Cnt_switch==1){
      
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
       
      /* CdOrbs_Grid */
 
      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){

        al = -1;
	for (L0=0; L0<=Spe_MaxL_Basis[Cwan]; L0++){
	  for (Mul0=0; Mul0<Spe_Num_CBasis[Cwan][L0]; Mul0++){
	    for (M0=0; M0<=2*L0; M0++){

	      al++;
	      tmp0 = 0.0;
	      tmpx = 0.0;
	      tmpy = 0.0;
	      tmpz = 0.0;

	      for (p=0; p<Spe_Specified_Num[Cwan][al]; p++){
	        j = Spe_Trans_Orbital[Cwan][al][p];  
	        tmp0 += CntCoes[Mc_AN][al][p]*Orbs_Grid[Mc_AN][Nc][j];/* AITUNE */
	        tmpx += CntCoes[Mc_AN][al][p]*dOrbs_Grid[0][Mc_AN][Nc][j];/* AITUNE */
	        tmpy += CntCoes[Mc_AN][al][p]*dOrbs_Grid[1][Mc_AN][Nc][j];/* AITUNE */
	        tmpz += CntCoes[Mc_AN][al][p]*dOrbs_Grid[2][Mc_AN][Nc][j];/* AITUNE */
	      }

	      COrbs_Grid[Mc_AN][al][Nc] = (Type_Orbs_Grid)tmp0;
	      CdOrbs_Grid[0][Mc_AN][al][Nc] = (Type_Orbs_Grid)tmpx;
	      CdOrbs_Grid[1][Mc_AN][al][Nc] = (Type_Orbs_Grid)tmpy;
	      CdOrbs_Grid[2][Mc_AN][al][Nc] = (Type_Orbs_Grid)tmpz;
	    }
	  }
        }
      }

    }

    /**********************************************
     MPI:

     CdOrbs_Grid    
    ***********************************************/

    /* allocation of arrays  */
    Snd_Size = (int*)malloc(sizeof(int)*numprocs); 
    Rcv_Size = (int*)malloc(sizeof(int)*numprocs); 

    /* find data size for sending and receiving */

    My_Max = -10000;
    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){
        /*  sending size */
        if (F_Snd_Num[IDS]!=0){
          /* find data size */ 
          size1 = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            size1 += GridN_Atom[Gc_AN]*Spe_Total_CNO[Cwan];
          }

          Snd_Size[IDS] = size1;
          MPI_Isend(&size1, 1, MPI_INT, IDS, tag, mpi_comm_level1, &request);
        }
        else{
          Snd_Size[IDS] = 0;
        }

        /*  receiving size */
        if (F_Rcv_Num[IDR]!=0){
          MPI_Recv(&size2, 1, MPI_INT, IDR, tag, mpi_comm_level1, &stat);
          Rcv_Size[IDR] = size2;
        }
        else{
          Rcv_Size[IDR] = 0;
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
      else{
        Snd_Size[IDS] = 0;
        Rcv_Size[IDR] = 0;
      }

      if (My_Max<Snd_Size[IDS]) My_Max = Snd_Size[IDS];
      if (My_Max<Rcv_Size[IDR]) My_Max = Rcv_Size[IDR];

    }  

    MPI_Allreduce(&My_Max, &Max_Size, 1, MPI_INT, MPI_MAX, mpi_comm_level1);
    /* allocation of arrays */ 
    tmp_array    = (double*)malloc(sizeof(double)*Max_Size);
    tmpx_array   = (double*)malloc(sizeof(double)*Max_Size);
    tmpy_array   = (double*)malloc(sizeof(double)*Max_Size);
    tmpz_array   = (double*)malloc(sizeof(double)*Max_Size);
    tmp_array2   = (double*)malloc(sizeof(double)*Max_Size);
    tmpx_array2  = (double*)malloc(sizeof(double)*Max_Size);
    tmpy_array2  = (double*)malloc(sizeof(double)*Max_Size);
    tmpz_array2  = (double*)malloc(sizeof(double)*Max_Size);

    /* send and receive COrbs_Grid */

    for (ID=0; ID<numprocs; ID++){

      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;

      if (ID!=0){

        /* sending of data */ 

        if (F_Snd_Num[IDS]!=0){

          /* find data size */
          size1 = Snd_Size[IDS];

          /* multidimentional array to vector array */
          k = 0; 
          for (n=0; n<F_Snd_Num[IDS]; n++){
            Mc_AN = Snd_MAN[IDS][n];
            Gc_AN = Snd_GAN[IDS][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 
            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                tmp_array[k] = COrbs_Grid[Mc_AN][i][Nc];
                tmpx_array[k] = CdOrbs_Grid[0][Mc_AN][i][Nc];
                tmpy_array[k] = CdOrbs_Grid[1][Mc_AN][i][Nc];
                tmpz_array[k] = CdOrbs_Grid[2][Mc_AN][i][Nc];
                k++;
              }          
            }
          } 

          /* MPI_Isend */
          MPI_Isend(&tmp_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
          MPI_Isend(&tmpx_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
          MPI_Isend(&tmpy_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
          MPI_Isend(&tmpz_array[0], size1, MPI_DOUBLE, IDS,
                    tag, mpi_comm_level1, &request);
        }

        /* receiving of block data */

        if (F_Rcv_Num[IDR]!=0){

          /* find data size */
          size2 = Rcv_Size[IDR]; 

          /* MPI_Recv */
          MPI_Recv(&tmp_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);
          MPI_Recv(&tmpx_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);
          MPI_Recv(&tmpy_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);
          MPI_Recv(&tmpz_array2[0], size2, MPI_DOUBLE, IDR, tag, mpi_comm_level1, &stat);

          k = 0;
          Mc_AN = F_TopMAN[IDR] - 1;
          for (n=0; n<F_Rcv_Num[IDR]; n++){
            Mc_AN++;
            Gc_AN = Rcv_GAN[IDR][n];
            Cwan = WhatSpecies[Gc_AN];
            NO0 = Spe_Total_CNO[Cwan]; 

            for (i=0; i<NO0; i++){
              for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
                COrbs_Grid[Mc_AN][i][Nc] = tmp_array2[k];
                CdOrbs_Grid[0][Mc_AN][i][Nc] = tmpx_array2[k];
                CdOrbs_Grid[1][Mc_AN][i][Nc] = tmpy_array2[k];
                CdOrbs_Grid[2][Mc_AN][i][Nc] = tmpz_array2[k];
                k++;
              }          
            }
          }
        }
        if (F_Snd_Num[IDS]!=0) MPI_Wait(&request,&stat);
      } 
    }  

    /* freeing of arrays  */
    free(tmp_array);
    free(tmpx_array);
    free(tmpy_array);
    free(tmpz_array);
    free(tmp_array2);
    free(tmpx_array2);
    free(tmpy_array2);
    free(tmpz_array2);
    free(Snd_Size);
    free(Rcv_Size);
  }

  /**********************************************
              calculate Tmp_dDen_Grid
  ***********************************************/
    
  
  /* AITUNE ========================== */ 
  int OneD_Nloop = 0;
  int ai_MaxNc = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    int Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
      if(ai_MaxNc < GridN_Atom[Gc_AN]) {ai_MaxNc = GridN_Atom[Gc_AN];}
    }
  }  
  /* ai_MaxNc is maximum of GridN_Atom[] */

  int gNthrds;
#pragma omp parallel
  {
    gNthrds = omp_get_num_threads();
  }

  double**** ai_tmpDG_all = (double****)malloc(sizeof(double***)*gNthrds);
	
  /* ========================== AITUNE */ 

#pragma omp parallel shared(myid,G2ID,Orbs_Grid_FNAN,dOrbs_Grid_FNAN,List_YOUSO,Tmp_dDen_Grid,Orbs_Grid,dOrbs_Grid,COrbs_Grid,CdOrbs_Grid,Cnt_switch,Cnt_kind,GListTAtoms2,GListTAtoms1,NumOLG,CDM,SpinP_Grid_size,SpinP_switch,WhatSpecies,ncn,F_G2M,natn,Spe_Total_CNO,M2G) private(OMPID,Nthrds,Mc_AN,h_AN,Stime_atom,Etime_atom,Gc_AN,Cwan,NO0,Gh_AN,Mh_AN,Rnh,Hwan,NO1,spin,i,j,k,tmp_CDM,Nog,Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3,orbs0_0,orbs0_1,orbs0_2,orbs0_3,orbs1_0,orbs1_1,orbs1_2,orbs1_3,dorbs0_0,dorbs0_1,dorbs0_2,dorbs0_3,dorbs1_0,dorbs1_1,dorbs1_2,dorbs1_3,sum_0,sum_1,sum_2,sum_3,sum2_0,sum2_1,sum2_2,sum2_3,tmp0_0,tmp0_1,tmp0_2,tmp0_3,Nc,Nh,orbs0,orbs1,dorbs0,dorbs1,sum,sum2,tmp0)
  {

//    orbs0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    dorbs0 = (double**)malloc(sizeof(double*)*3);
    for (k=0; k<3; k++){
      dorbs0[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
    dorbs1 = (double**)malloc(sizeof(double*)*3);
    for (k=0; k<3; k++){
      dorbs1[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }

//    orbs0_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs0_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs0_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs0_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs1_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs1_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs1_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
//    orbs1_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    dorbs0_0 = (double**)malloc(sizeof(double*)*3);
    dorbs0_1 = (double**)malloc(sizeof(double*)*3);
    dorbs0_2 = (double**)malloc(sizeof(double*)*3);
    dorbs0_3 = (double**)malloc(sizeof(double*)*3);
    for (k=0; k<3; k++){
      dorbs0_0[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs0_1[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs0_2[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs0_3[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }
    dorbs1_0 = (double**)malloc(sizeof(double*)*3);
    dorbs1_1 = (double**)malloc(sizeof(double*)*3);
    dorbs1_2 = (double**)malloc(sizeof(double*)*3);
    dorbs1_3 = (double**)malloc(sizeof(double*)*3);
    for (k=0; k<3; k++){
      dorbs1_0[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs1_1[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs1_2[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
      dorbs1_3[k] = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    }

    tmp_CDM = (double***)malloc(sizeof(double**)*(SpinP_switch+1)); 
    for (i=0; i<(SpinP_switch+1); i++){
      tmp_CDM[i] = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
      for (j=0; j<List_YOUSO[7]; j++){
      	tmp_CDM[i][j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
      }
    }

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();

	
    /* AITUNE ========================== */  


    double **ai_tmpDGs[7];
    {
      for (k=0; k<7; k++){
         ai_tmpDGs[k] = (double**)malloc(sizeof(double*)*(SpinP_Grid_size+1));
         int spin;
         for (spin=0; spin<(SpinP_Grid_size+1); spin++){
         	ai_tmpDGs[k][spin] = (double*)malloc(sizeof(double)* ai_MaxNc);
         }
      }
    }
    ai_tmpDG_all[OMPID] = ai_tmpDGs;
    /* ==================================== AITUNE */


    /* for (Mc_AN=(OMPID+1); Mc_AN<=Matomnum; Mc_AN+=Nthrds) AITUNE */
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 
	  
      for (k=0; k<7; k++){
         int spin;
         for (spin=0; spin<(SpinP_Grid_size+1); spin++){
         	int Nc;
         	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
         	  ai_tmpDGs[k][spin][Nc] = 0.0;
         	}
         }
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	/* set data on h_AN */
    
	Gh_AN = natn[Gc_AN][h_AN];
	Mh_AN = F_G2M[Gh_AN];
	Rnh = ncn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	NO1 = Spe_Total_CNO[Hwan];

	/* store CDM into tmp_CDM */

	for (spin=0; spin<SpinP_switch+1; spin++){
	  for (i=0; i<NO0; i++){
	    for (j=0; j<NO1; j++){
	      tmp_CDM[spin][i][j] = CDM[spin][Mc_AN][h_AN][i][j];
	    }
	  }
	}

	/* summation of non-zero elements */
	/* for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++) */
#pragma omp for
	for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]-3; Nog+=4){

	  Nc_0 = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nc_1 = GListTAtoms1[Mc_AN][h_AN][Nog+1];
	  Nc_2 = GListTAtoms1[Mc_AN][h_AN][Nog+2];
	  Nc_3 = GListTAtoms1[Mc_AN][h_AN][Nog+3];
	  
	  Nh_0 = GListTAtoms2[Mc_AN][h_AN][Nog];
	  Nh_1 = GListTAtoms2[Mc_AN][h_AN][Nog+1];
	  Nh_2 = GListTAtoms2[Mc_AN][h_AN][Nog+2];
	  Nh_3 = GListTAtoms2[Mc_AN][h_AN][Nog+3];
	  
	  /* Now under the orbital optimization */
	  if (Cnt_kind==0 && Cnt_switch==1){
	    //for (i=0; i<NO0; i++){
	    //  orbs0_0[i] = COrbs_Grid[Mc_AN][i][Nc_0];
	    //  orbs0_1[i] = COrbs_Grid[Mc_AN][i][Nc_1];
	    //  orbs0_2[i] = COrbs_Grid[Mc_AN][i][Nc_2];
	    //  orbs0_3[i] = COrbs_Grid[Mc_AN][i][Nc_3];
	    //}
	    for (k=0; k<3; k++){
	      for (i=0; i<NO0; i++){
	         dorbs0_0[k][i] = CdOrbs_Grid[k][Mc_AN][i][Nc_0];
	         dorbs0_1[k][i] = CdOrbs_Grid[k][Mc_AN][i][Nc_1];
	         dorbs0_2[k][i] = CdOrbs_Grid[k][Mc_AN][i][Nc_2];
	         dorbs0_3[k][i] = CdOrbs_Grid[k][Mc_AN][i][Nc_3];
	      }
	    }
	    //for (j=0; j<NO1; j++){
	    //  orbs1_0[j] = COrbs_Grid[Mh_AN][j][Nh_0];
	    //  orbs1_1[j] = COrbs_Grid[Mh_AN][j][Nh_1];
	    //  orbs1_2[j] = COrbs_Grid[Mh_AN][j][Nh_2];
	    //  orbs1_3[j] = COrbs_Grid[Mh_AN][j][Nh_3];
	    //}
	    for (k=0; k<3; k++){
	      for (j=0; j<NO1; j++){
	         dorbs1_0[k][j] = CdOrbs_Grid[k][Mh_AN][j][Nh_0];
	         dorbs1_1[k][j] = CdOrbs_Grid[k][Mh_AN][j][Nh_1];
	         dorbs1_2[k][j] = CdOrbs_Grid[k][Mh_AN][j][Nh_2];
	         dorbs1_3[k][j] = CdOrbs_Grid[k][Mh_AN][j][Nh_3];
	      }
	    }
	  }
	  /* else if ! "now under the orbital optimization" */
	  else{
	    //for (i=0; i<NO0; i++){
	    //  orbs0_0[i] = Orbs_Grid[Mc_AN][Nc_0][i];
	    //  orbs0_1[i] = Orbs_Grid[Mc_AN][Nc_1][i];
	    //  orbs0_2[i] = Orbs_Grid[Mc_AN][Nc_2][i];
	    //  orbs0_3[i] = Orbs_Grid[Mc_AN][Nc_3][i]; 
       //}
	    for (k=0; k<3; k++){
	      for (i=0; i<NO0; i++){
		      dorbs0_0[k][i] = dOrbs_Grid[k][Mc_AN][Nc_0][i];
		      dorbs0_1[k][i] = dOrbs_Grid[k][Mc_AN][Nc_1][i];
		      dorbs0_2[k][i] = dOrbs_Grid[k][Mc_AN][Nc_2][i];
		      dorbs0_3[k][i] = dOrbs_Grid[k][Mc_AN][Nc_3][i]; 
	      }
	    }

       if (G2ID[Gh_AN]==myid){
	      //for (j=0; j<NO1; j++){
	      //   orbs1_0[j] = Orbs_Grid[Mh_AN][Nh_0][j];
	      //   orbs1_1[j] = Orbs_Grid[Mh_AN][Nh_1][j];
	      //   orbs1_2[j] = Orbs_Grid[Mh_AN][Nh_2][j];
	      //   orbs1_3[j] = Orbs_Grid[Mh_AN][Nh_3][j]; 
         //}
	      for (k=0; k<3; k++){
	         for (j=0; j<NO1; j++){
		         dorbs1_0[k][j] = dOrbs_Grid[k][Mh_AN][Nh_0][j];
		         dorbs1_1[k][j] = dOrbs_Grid[k][Mh_AN][Nh_1][j];
		         dorbs1_2[k][j] = dOrbs_Grid[k][Mh_AN][Nh_2][j];
		         dorbs1_3[k][j] = dOrbs_Grid[k][Mh_AN][Nh_3][j]; 
	         }
	      }
	    }
       else{
	      //for (j=0; j<NO1; j++){
	      //   orbs1_0[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog  ][j];
	      //   orbs1_1[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+1][j];
	      //   orbs1_2[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+2][j];
	      //   orbs1_3[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+3][j]; 
         //}
	      for (k=0; k<3; k++){
	         for (j=0; j<NO1; j++){
		         dorbs1_0[k][j] = dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nog  ][j];
		         dorbs1_1[k][j] = dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nog+1][j];
		         dorbs1_2[k][j] = dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nog+2][j];
		         dorbs1_3[k][j] = dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nog+3][j]; 
	         }
	      }
	    }
	  }
	  

	    /* Tmp_dDen_Grid */

     for (k=0; k<3; k++){
       int itmp,itmp1,itmp2;
       itmp = (k+2)%3;
       itmp1 = k+1;
       itmp2 = itmp+4;
	    for (spin=0; spin<=SpinP_switch; spin++){

	      sum_0 = 0.0;  sum2_0 = 0.0;
	      sum_1 = 0.0;  sum2_1 = 0.0;
	      sum_2 = 0.0;  sum2_2 = 0.0;
	      sum_3 = 0.0;  sum2_3 = 0.0;

	      for (i=0; i<NO0; i++){
	          tmp0_0 = 0.0;
	          tmp0_1 = 0.0;
	          tmp0_2 = 0.0;
	          tmp0_3 = 0.0;

	        for (j=0; j<NO1; j++){
		        tmp0_0 += dorbs1_0[k][j]*tmp_CDM[spin][i][j];
		        tmp0_1 += dorbs1_1[k][j]*tmp_CDM[spin][i][j];
		        tmp0_2 += dorbs1_2[k][j]*tmp_CDM[spin][i][j];
		        tmp0_3 += dorbs1_3[k][j]*tmp_CDM[spin][i][j];
	        }

	           sum_0 += dorbs0_0[k][i]*tmp0_0;
	           sum_1 += dorbs0_1[k][i]*tmp0_1;
	           sum_2 += dorbs0_2[k][i]*tmp0_2;
	           sum_3 += dorbs0_3[k][i]*tmp0_3;
	           sum2_0 += dorbs0_0[itmp][i]*tmp0_0;
	           sum2_1 += dorbs0_1[itmp][i]*tmp0_1;
	           sum2_2 += dorbs0_2[itmp][i]*tmp0_2;
	           sum2_3 += dorbs0_3[itmp][i]*tmp0_3;
	      }
		
	      ai_tmpDGs[itmp1][spin][Nc_0] += sum_0;
	      ai_tmpDGs[itmp1][spin][Nc_1] += sum_1;
	      ai_tmpDGs[itmp1][spin][Nc_2] += sum_2;
	      ai_tmpDGs[itmp1][spin][Nc_3] += sum_3;
	      ai_tmpDGs[itmp2][spin][Nc_0] += sum2_0;
	      ai_tmpDGs[itmp2][spin][Nc_1] += sum2_1;
	      ai_tmpDGs[itmp2][spin][Nc_2] += sum2_2;
	      ai_tmpDGs[itmp2][spin][Nc_3] += sum2_3;
	      /*
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_0] += sum_0;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_1] += sum_1;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_2] += sum_2;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_3] += sum_3;
	      */
       } /* spin */
	  } /* k */
     /* kokomade 2.May.2015 */

     for (k=0; k<3; k++){
       int itmp,itmp1,itmp2;
       itmp = (k+2)%3;
       itmp1 = k+1;
       itmp2 = itmp+4;
	    for (spin=2; spin<=SpinP_switch; spin++){

	      sum_0 = 0.0;  sum2_0 = 0.0;
	      sum_1 = 0.0;  sum2_1 = 0.0;
	      sum_2 = 0.0;  sum2_2 = 0.0;
	      sum_3 = 0.0;  sum2_3 = 0.0;

	      for (i=0; i<NO0; i++){

	          tmp0_0 = 0.0;
	          tmp0_1 = 0.0;
	          tmp0_2 = 0.0;
	          tmp0_3 = 0.0;

	        for (j=0; j<NO1; j++){
		        tmp0_0 += dorbs1_0[k][j]*tmp_CDM[spin][i][j];
		        tmp0_1 += dorbs1_1[k][j]*tmp_CDM[spin][i][j];
		        tmp0_2 += dorbs1_2[k][j]*tmp_CDM[spin][i][j];
		        tmp0_3 += dorbs1_3[k][j]*tmp_CDM[spin][i][j];
	        }

	        sum_0 += dorbs0_0[k][i]*tmp0_0;
	        sum_1 += dorbs0_1[k][i]*tmp0_1;
	        sum_2 += dorbs0_2[k][i]*tmp0_2;
	        sum_3 += dorbs0_3[k][i]*tmp0_3;
	        sum2_0 += dorbs0_0[itmp][i]*tmp0_0;
	        sum2_1 += dorbs0_1[itmp][i]*tmp0_1;
	        sum2_2 += dorbs0_2[itmp][i]*tmp0_2;
	        sum2_3 += dorbs0_3[itmp][i]*tmp0_3;
	      }
	      ai_tmpDGs[itmp1][spin+2][Nc_0] += sum_0;
	      ai_tmpDGs[itmp1][spin+2][Nc_1] += sum_1;
	      ai_tmpDGs[itmp1][spin+2][Nc_2] += sum_2;
	      ai_tmpDGs[itmp1][spin+2][Nc_3] += sum_3;
	      ai_tmpDGs[itmp2][spin][Nc_0] += sum2_0;
	      ai_tmpDGs[itmp2][spin][Nc_1] += sum2_1;
	      ai_tmpDGs[itmp2][spin][Nc_2] += sum2_2;
	      ai_tmpDGs[itmp2][spin][Nc_3] += sum2_3;
	      /*
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_0] += sum_0;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_1] += sum_1;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_2] += sum_2;
	        Tmp_dDen_Grid[k][spin][Mc_AN][Nc_3] += sum_3;
	      */
       } /* spin */
	  } /* k */

	} /* Nog */

#pragma omp for
	for (Nog = NumOLG[Mc_AN][h_AN] - (NumOLG[Mc_AN][h_AN] % 4); Nog<NumOLG[Mc_AN][h_AN]; Nog++){
	  /*for (; Nog<NumOLG[Mc_AN][h_AN]; Nog++)*/
	
	  Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
 

	  if (Cnt_kind==0 && Cnt_switch==1){
	    //for (i=0; i<NO0; i++){
	    //  orbs0[i] = COrbs_Grid[Mc_AN][i][Nc];
	    //}
	    for (k=0; k<3; k++){
	      for (i=0; i<NO0; i++){
	         dorbs0[k][i] = CdOrbs_Grid[k][Mc_AN][i][Nc];
	      }
	    }
	    //for (j=0; j<NO1; j++){
	    //  orbs1[j] = COrbs_Grid[Mh_AN][j][Nh];
	    //}
	    for (k=0; k<3; k++){
	      for (j=0; j<NO1; j++){
	         dorbs1[k][j] = CdOrbs_Grid[k][Mh_AN][j][Nh];
	      }
	    }
	  }
	  else{
	    //for (i=0; i<NO0; i++){
	    //  orbs0[i] = Orbs_Grid[Mc_AN][Nc][i];
	    //}
	    for (k=0; k<3; k++){
	      for (i=0; i<NO0; i++){
	         dorbs0[k][i] = dOrbs_Grid[k][Mc_AN][Nc][i];
	      }
	    }

	    if (G2ID[Gh_AN]==myid){
	      //for (j=0; j<NO1; j++){
	      //  orbs1[j] = Orbs_Grid[Mh_AN][j][Nh];
	      //}
	      for (k=0; k<3; k++){
	         for (j=0; j<NO1; j++){
		         dorbs1[k][j] = dOrbs_Grid[k][Mh_AN][Nh][j];
	         }
	      }
	    }
	    else{
	      //for (j=0; j<NO1; j++){
	      //  orbs1[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];
	      //}
	      for (k=0; k<3; k++){
	         for (j=0; j<NO1; j++){
		         dorbs1[k][j] = dOrbs_Grid_FNAN[k][Mc_AN][h_AN][Nog][j];
	         }
	      }
	    }
	  }

	  for (spin=0; spin<=SpinP_switch; spin++){
 
	    for (k=0; k<3; k++){
         int itmp,itmp1,itmp2;
         itmp = (k+2)%3;
         itmp1 = k+1;
         itmp2 = itmp+4;
	      sum = 0.0; sum2 = 0.0;
	      for (i=0; i<NO0; i++){
	        tmp0 = 0.0;
	        for (j=0; j<NO1; j++){
           	tmp0 += dorbs1[k][j]*tmp_CDM[spin][i][j];
	        }
	        sum += dorbs0[k][i]*tmp0;
	        sum2 += dorbs0[itmp][i]*tmp0;
	      }
 
	      ai_tmpDGs[itmp1][spin][Nc] += sum;
	      ai_tmpDGs[itmp2][spin][Nc] += sum2;
	      /*Tmp_dDen_Grid[k][spin][Mc_AN][Nc] += sum;*/
       }
	  } /* spin */

	  for (spin=2; spin<=SpinP_switch; spin++){
 
	    for (k=0; k<3; k++){
         int itmp,itmp1,itmp2;
         itmp = (k+2)%3;
         itmp1 = k+1;
         itmp2 = itmp+4;
	      sum = 0.0; sum2 = 0.0;
	      for (i=0; i<NO0; i++){
	        tmp0 = 0.0;
	        for (j=0; j<NO1; j++){
           	tmp0 += dorbs1[k][j]*tmp_CDM[spin][i][j];
	        }
	        sum += dorbs0[k][i]*tmp0;
	        sum2 += dorbs0[itmp][i]*tmp0;
	      }
 
	      ai_tmpDGs[itmp1][spin+2][Nc] += sum;
	      ai_tmpDGs[itmp2][spin+2][Nc] += sum2;
	      /*Tmp_dDen_Grid[k][spin][Mc_AN][Nc] += sum;*/
       }
	  } /* spin */

	} /* Nog */
	
      } /* h_AN */

      /* AITUNE   merge temporary buffer for all omp threads */	
	   for (k=0; k<7; k++){
         for (spin=0; spin<=(SpinP_Grid_size); spin++){
	      int Nc;
#pragma omp for
         	for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
         	  double sum = 0.0;
         	  int th;
         	  for(th = 0; th < Nthrds; th++){
         	    sum += ai_tmpDG_all[th][k][spin][Nc];
         	  }
         	  Tmp_dDen_Grid[k][spin][Mc_AN][Nc] += sum;
         	}
         }
		}

    } /* Mc_AN */

    /* freeing of arrays */ 

    //free(orbs0);
    //free(orbs1);
	 for (k=0; k<3; k++){
      free(dorbs0[k]);
      free(dorbs1[k]);
    }
    free(dorbs0);
    free(dorbs1);

    //free(orbs0_0);
    //free(orbs0_1);
    //free(orbs0_2);
    //free(orbs0_3);
    //free(orbs1_0);
    //free(orbs1_1);
    //free(orbs1_2);
    //free(orbs1_3);
	 for (k=0; k<3; k++){
      free(dorbs0_0[k]);
      free(dorbs0_1[k]);
      free(dorbs0_2[k]);
      free(dorbs0_3[k]);
      free(dorbs1_0[k]);
      free(dorbs1_1[k]);
      free(dorbs1_2[k]);
      free(dorbs1_3[k]);
    }
    free(dorbs0_0);
    free(dorbs0_1);
    free(dorbs0_2);
    free(dorbs0_3);
    free(dorbs1_0);
    free(dorbs1_1);
    free(dorbs1_2);
    free(dorbs1_3);

    for (i=0; i<(SpinP_switch+1); i++){
      for (j=0; j<List_YOUSO[7]; j++){
      	free(tmp_CDM[i][j]);
      }
      free(tmp_CDM[i]);
    }
    free(tmp_CDM);

	 for (k=0; k<7; k++){
      for (i=0; i<(SpinP_Grid_size+1); i++){
         free(ai_tmpDGs[k][i]); /* AITUNE */
      }
      free(ai_tmpDGs[k]);
    }

#pragma omp flush(Tmp_dDen_Grid)

  } /* #pragma omp parallel */
  
  free(ai_tmpDG_all);

  /******************************************************
      MPI communication from the partitions A to B 
  ******************************************************/
  
  /* copy Tmp_dDen_Grid to Den_Snd_Grid_A2B */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;
  
  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

	   for (k=0; k<7; k++){
         int itmp;
         itmp = k*(SpinP_Grid_size+1);
         for (i=0; i<(SpinP_Grid_size+1); i++){
            Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*(SpinP_Grid_size+1)*7+itmp+i] = Tmp_dDen_Grid[k][i][Mc_AN][AN];
         } /* i SponP */
      } /* k */
      /*
      if (SpinP_switch==0){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]] = Tmp_Den_Grid[0][Mc_AN][AN];
      }
      else if (SpinP_switch==1){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*2+1] = Tmp_Den_Grid[1][Mc_AN][AN];
      }
      else if (SpinP_switch==3){
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+0] = Tmp_Den_Grid[0][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+1] = Tmp_Den_Grid[1][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+2] = Tmp_Den_Grid[2][Mc_AN][AN];
        Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]*4+3] = Tmp_Den_Grid[3][Mc_AN][AN];
      }
      */

      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: A to B */  

  request_send = malloc(sizeof(MPI_Request)*NN_A2B_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_A2B_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_A2B_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_A2B_R);

  NN_S = 0;
  NN_R = 0;

  tag = 999;
  for (ID=1; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Den_Snd_Grid_A2B[IDS][0], Num_Snd_Grid_A2B[IDS]*(SpinP_Grid_size+1)*7, 
	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Irecv( &Den_Rcv_Grid_A2B[IDR][0], Num_Rcv_Grid_A2B[IDR]*(SpinP_Grid_size+1)*7, 
  	         MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* for myid */
  for (i=0; i<Num_Rcv_Grid_A2B[myid]*(SpinP_Grid_size+1)*7; i++){
    Den_Rcv_Grid_A2B[myid][i] = Den_Snd_Grid_A2B[myid][i];
  }

  /******************************************************
   superposition of rho_i to calculate charge density 
   in the partition B.
  ******************************************************/

  /* initialize arrays */

   for (k=0; k<7; k++){
      for (spin=0; spin<(SpinP_Grid_size+1); spin++){
        for (BN=0; BN<My_NumGridB_AB; BN++){
          dd2Density_Grid_B[k][spin][BN] = 0.0;
        }
      }
   }

  /* superposition of densities rho_i */

  for (ID=0; ID<numprocs; ID++){

    for (LN=0; LN<Num_Rcv_Grid_A2B[ID]; LN++){

      BN    = Index_Rcv_Grid_A2B[ID][3*LN+0];      
      Gc_AN = Index_Rcv_Grid_A2B[ID][3*LN+1];        
      GRc   = Index_Rcv_Grid_A2B[ID][3*LN+2]; 

      if (Solver!=4 || (Solver==4 && atv_ijk[GRc][1]==0 )){

	      for (k=0; k<7; k++){
            int itmp;
            itmp = k*(SpinP_Grid_size+1);
            for (i=0; i<(SpinP_Grid_size+1); i++){
	            dd2Density_Grid_B[k][i][BN] += Den_Rcv_Grid_A2B[ID][LN*(SpinP_Grid_size+1)*7+itmp+i];
            }
         }

      } /* if (Solve!=4.....) */           

    } /* AN */ 
  } /* ID */  


  /****************************************************
   Conjugate complex of Density_Grid[3][MN] due to
   difference in the definition between density matrix
   and charge density
  ****************************************************/

//  for(k=0; k<3; k++){
//     for (BN=0; BN<My_NumGridB_AB; BN++){
//       dDensity_Grid_B[k][3][BN] = -dDensity_Grid_B[k][3][BN]; 
//       dDensity_Grid_B[k][5][BN] = -dDensity_Grid_B[k][5][BN]; 
//     }
//  }


/*******************************************************
 double ***dd2Density_Grid_B; 
  electron densities on grids in the partition B
  size: dd2Density_Grid[7][6][My_NumGridB_AB]
       (n_dd2)_{ab}^0 = 0
       (n_dd2)_{ab}^1 = (\partial_x \psi)^{a*} \partial_x \psi^b
       (n_dd2)_{ab}^2 = (\partial_y \psi)^{a*} \partial_y \psi^b
       (n_dd2)_{ab}^3 = (\partial_z \psi)^{a*} \partial_z \psi^b
       (n_dd2)_{ab}^4 = (\partial_x \psi)^{a*} \partial_y \psi^b
       (n_dd2)_{ab}^5 = (\partial_y \psi)^{a*} \partial_z \psi^b
       (n_dd2)_{ab}^6 = (\partial_z \psi)^{a*} \partial_x \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       dd2Density_Grid_B[k][0][BN]　:  Re[nd_{11}^k]
       dd2Density_Grid_B[k][1][BN]　:  Re[nd_{22}^k] 
       dd2Density_Grid_B[k][2][BN]　:  Re[nd_{12}^k]
       dd2Density_Grid_B[k][3][BN]　:  Im[nd_{12}^k]
       dd2Density_Grid_B[k][4][BN]　:  Re[nd_{21}^k]
       dd2Density_Grid_B[k][5][BN]　: -Im[nd_{21}^k]
*******************************************************/


//  /******************************************************
//             MPI: from the partitions B to D
//  ******************************************************/
//
//  Density_Grid_Copy_B2D();

  /* freeing of arrays */

  for (k=0; k<7; k++){
     for (i=0; i<(SpinP_Grid_size+1); i++){
       for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
         free(Tmp_dDen_Grid[k][i][Mc_AN]);
       }
       free(Tmp_dDen_Grid[k][i]);
     }
     free(Tmp_dDen_Grid[k]);
  }
  free(Tmp_dDen_Grid);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Snd_Grid_A2B[ID]);
  }  
  free(Den_Snd_Grid_A2B);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Rcv_Grid_A2B[ID]);
  }
  free(Den_Rcv_Grid_A2B);

  return time0;
}

