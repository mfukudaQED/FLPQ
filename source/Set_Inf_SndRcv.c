/**********************************************************************
  Set_inf_SndRcv.c:

     Set_inf_SndRcv.c is a subrutine to construct MPI structures.

  Log of Set_inf_SndRcv.c:

     10/May/2019  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "omp.h"
#include "flpq.h"


void Set_Inf_SndRcv()
{ 
  int i,ID,IDS,IDR,Mc_AN,Gc_AN,Num,ID1,Lh_AN,Gh_AN;
  int myid,numprocs,tag=999;
  int *flag_DoubleCounting;
  int **Rcv_FGAN;
  int **Snd_FGAN;

  int Rn,Rn2,m1,m2,m3,n1,n2,n3,j,po,Gj_AN;
  MPI_Status stat;
  MPI_Request request;

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /*********************************
   allocation of arrays:

   int flag_DoubleCounting[atomnum+1]
  *********************************/

  flag_DoubleCounting = (int*)malloc(sizeof(int)*(atomnum+1));

  /*********************************
   initialize

      F_Snd_Num
      F_Rcv_Num
  *********************************/

  for (ID=0; ID<numprocs; ID++){
    F_Snd_Num[ID] = 0;
    F_Rcv_Num[ID] = 0;
  }
    
  /************************************************
      find F_Rcv_Num and S_Rcv_Num
  *************************************************/

  /* initialize flag_DoubleCounting */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    flag_DoubleCounting[Gc_AN] = 0;
  }

  //for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
  //        printf("Gc_AN, ID = %d, %d\n", Gc_AN, G2ID[Gc_AN]);
  //}

  /* find F_Rcv_Num */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      //printf("%d,%d,%d,%d\n",Matomnum,Mc_AN,FNAN[Gc_AN],Lh_AN);
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){
        F_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
      }
    }
  }

  /************************************************
   allocation of array:

   int Rcv_GAN[numprocs]
              [F_Rcv_Num[ID]]
  *************************************************/

  Rcv_GAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_GAN[ID] = (int*)malloc(sizeof(int)*(F_Rcv_Num[ID]));
  }

  Rcv_FGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Rcv_FGAN[ID] = (int*)malloc(sizeof(int)*F_Rcv_Num[ID]);
  }

  /************************************************
             set Rcv_FGAN and Rcv_SGAN
  *************************************************/
  
  /* initialize F_Rcv_Num */

  for (ID=0; ID<numprocs; ID++){
    F_Rcv_Num[ID] = 0;
  }

  /* initialized flag_DoubleCounting */

  for (Gc_AN=1; Gc_AN<=atomnum; Gc_AN++){
    flag_DoubleCounting[Gc_AN] = 0;
  }

  /* set Rcv_FGAN */

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    for (Lh_AN=0; Lh_AN<=FNAN[Gc_AN]; Lh_AN++){
      Gh_AN = natn[Gc_AN][Lh_AN];
      ID1 = G2ID[Gh_AN];
      if (flag_DoubleCounting[Gh_AN]==0 && ID1!=myid){

        Rcv_FGAN[ID1][F_Rcv_Num[ID1]] = Gh_AN;
        F_Rcv_Num[ID1]++;
        flag_DoubleCounting[Gh_AN] = 1;
      }
    }
  }

  /*****************************************
       MPI:  F_Rcv_Num
  *****************************************/

  /* MPI_Barrier */
  MPI_Barrier(mpi_comm_level1);
  
  tag = 999;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&F_Rcv_Num[IDS],1,MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &F_Snd_Num[IDR],1,MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  Snd_FGAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_FGAN[ID] = (int*)malloc(sizeof(int)*F_Snd_Num[ID]);
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=0){
      IDS = (myid + ID) % numprocs;
      IDR = (myid - ID + numprocs) % numprocs;
      MPI_Isend(&Rcv_FGAN[IDS][0],F_Rcv_Num[IDS],MPI_INT,IDS,tag,mpi_comm_level1,&request);
      MPI_Recv( &Snd_FGAN[IDR][0],F_Snd_Num[IDR],MPI_INT,IDR,tag,mpi_comm_level1,&stat); 
      MPI_Wait(&request,&stat);
    }
  }

  /********************************************************
    allocation of arrays:
  
    int Snd_MAN[numprocs][F_Snd_Num[ID]+1]
    int Snd_GAN[numprocs][F_Snd_Num[ID]+1]
  *********************************************************/

  Snd_MAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_MAN[ID] = (int*)malloc(sizeof(int)*(F_Snd_Num[ID]+1));
  }

  Snd_GAN = (int**)malloc(sizeof(int*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Snd_GAN[ID] = (int*)malloc(sizeof(int)*(F_Snd_Num[ID]+1));
  }

  /************************************************
      find data structures to send informations
      related to FNAN from myid to the other IDs 
  *************************************************/

  for (ID=0; ID<numprocs; ID++){

    Num = 0;

    for (i=0; i<F_Snd_Num[ID]; i++){ 

      for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

	Gc_AN = M2G[Mc_AN];

	if (Gc_AN==Snd_FGAN[ID][i]){

	  Snd_MAN[ID][Num] = Mc_AN;
	  Snd_GAN[ID][Num] = Gc_AN;
	  Num++;
	}      
      }
    }
  }

  /************************************************
   MPI:

     Snd_GAN
     Rcv_GAN
  *************************************************/

  for (ID=0; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (ID!=0){
      tag = 999;
      MPI_Isend(&Snd_GAN[IDS][0], F_Snd_Num[IDS],
                MPI_INT,IDS, tag, mpi_comm_level1, &request);
      MPI_Recv(&Rcv_GAN[IDR][0],  F_Rcv_Num[IDR],
                MPI_INT, IDR, tag, mpi_comm_level1, &stat);
      MPI_Wait(&request,&stat);
    }
  }

  /************************************************
          setting of F_TopMAN

    F_TopMAN give the first intermediate
    atom number in atoms sent from ID in the size of
    F_Rcv_Num[ID].
  *************************************************/

  Num = Matomnum + 1;
  for (ID=0; ID<numprocs; ID++){
    if (F_Rcv_Num[ID]!=0 && ID!=myid){
      F_TopMAN[ID] = Num;
      Num = Num + F_Rcv_Num[ID];
    }
  }


  /************************************************
       MatomnumF = the sum of F_Rcv_Num[ID]
  *************************************************/

  MatomnumF = 0;
  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid) MatomnumF += F_Rcv_Num[ID];
  }

  /************************************************
       allocation of arrays:
         
          F_M2G
          F_G2M
  *************************************************/

  F_M2G = (int*)malloc(sizeof(int)*(Matomnum+MatomnumF+1));
  F_G2M = (int*)malloc(sizeof(int)*(atomnum+1));

  /************************************************
           setting of F_G2M

    F_G2M and S_G2M give a conversion from the
    global atom number to the intermediate atom number
    for atoms sent from ID in the size of
    F_Rcv_Num[ID]. 
  *************************************************/
  
  /* initialization of F_G2M*/
  
  for (Gc_AN=0; Gc_AN<=atomnum; Gc_AN++){
    F_G2M[Gc_AN] = -1;
  } 
  
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = M2G[Mc_AN];
    F_G2M[Gc_AN] = Mc_AN;
    F_M2G[Mc_AN] = Gc_AN;
  }

  for (ID=0; ID<numprocs; ID++){
    if (ID!=myid && F_Rcv_Num[ID]!=0){
      for (Num=0; Num<F_Rcv_Num[ID]; Num++){
        Gc_AN = Rcv_GAN[ID][Num];
        F_G2M[Gc_AN] = F_TopMAN[ID] + Num;
        F_M2G[F_TopMAN[ID] + Num] = Gc_AN;  
      }
    }
  }

  /*********************************
   freeing of arrays:
  *********************************/

  for (ID=0; ID<numprocs; ID++){
    free(Snd_FGAN[ID]);
  }
  free(Snd_FGAN);

  for (ID=0; ID<numprocs; ID++){
    free(Rcv_FGAN[ID]);
  }
  free(Rcv_FGAN);

  free(flag_DoubleCounting);
} 

