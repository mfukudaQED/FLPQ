
/**********************************************************************
  main_flpq.c:

     main_flpq.c is the main routine of flpq.

  Log of main_flpq.c:

     27/Jun/2016  Released by M.Fukuda

***********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "mpi.h"
#include <omp.h>
#include "flpq.h"

int main(int argc, char *argv[]) 
{ 
  static int numprocs,myid;
  static int MD_iter,i,j,po,ip,n1;
  double TStime,TEtime;
  int kloop,mu;
  char filename[100];

  /* MPI initialize */

  mpi_comm_level1 = MPI_COMM_WORLD; 
  MPI_COMM_WORLD1 = MPI_COMM_WORLD; 

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  NUMPROCS_MPI_COMM_WORLD = numprocs;
  MYID_MPI_COMM_WORLD = myid;
  Num_Procs = numprocs;

  /* check argv */

  //if (argc==1){
  //  printf("\nCould not find an input file.\n\n");
  //  MPI_Finalize(); 
  //  exit(0);
  //} 


  if (myid==Host_ID){  
    printf("*******************************************************\n"); 
    printf("         Program code for LPQ calculations             \n");
    printf("                Written by M.FUKUDA                    \n");
    printf("*******************************************************\n\n"); 
  }

  /*****************************************************
               read input from flpq.inp
  *****************************************************/
  if (myid==Host_ID){  
         printf("start reading flpq.inp\n");
  }
  //read_input_flpq();
  sprintf(filename, "flpq.inp");
  Input_std(filename);

  //MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){  
    printf("end reading flpq.inp\n");
  }

  /*****************************************************
            read input from OpenMX_LPQ.dat
  *****************************************************/
  if (myid==Host_ID){  
    printf("start reading OpenMX_LPQ.bin\n");
  }

  read_input_openmx_lpq();

  //MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){  
    printf("end reading OpenMX_LPQ.bin\n");
    printf("\n");fflush(stdout);
  }

  //if ((Solver==4)&&(SpinP_switch<=1)){
  //  printf("Sorry. Non-relativistic NEGF is not available in flpq.\n");
  //  exit(0);
  //}
  //printf("OK");
  //exit(0);
    

  /*****************************************************
            read DM from OpenMX_LPQ_DM.bin
  *****************************************************/
    //if((flag_DMmu==0)&&(flag_DMkmu==0)){
    //  kloop=0;
    //  mu=0;
    //}
    //kloop=Num_DM_k;
    //mu=Num_DM_orb;
    
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){  
    printf("start read_input_DM\n");
  }

  read_input_DM( Num_DM_k, Num_DM_orb );

  if (myid==Host_ID){  
    printf("end read_input_DM\n");
    printf("\n");fflush(stdout);
  }


  /*****************************************************
        Set Grid_Origin_N to make a loop for xmesh 
  *****************************************************/
  Set_Grid_Origin_N();
  //printf("Ngrid1,Ngrid2,Ngrid3=%d, %d, %d\n",Ngrid1,Ngrid2,Ngrid3);fflush(stdout);
  //printf("flag_save_memory_3d=%d\n",flag_save_memory_3d);fflush(stdout);

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){  
    printf("Ngrid = %d, %d, %d\n",Ngrid1,Ngrid2,Ngrid3);
    for (i=1; i<=3; i++){
      printf("original lattice vector %d : %lf %lf %lf\n",i,tv_ori[i][1],tv_ori[i][2],tv_ori[i][3]);
    }
    for (i=1; i<=3; i++){
      printf("grid_vec %d : %lf %lf %lf\n",i,grid_vec[i][1],grid_vec[i][2],grid_vec[i][3]);
    }
    printf("start Set_Grid and out_LPQ\n");
    printf("\n");fflush(stdout);
  }

  if (flag_save_memory_3d==1){
    loop_Ngrid1 = Ngrid1;
    Ngrid1 = 1;
  }
  else{
    loop_Ngrid1 = 1;
  }

  /*****************************************************
        loop for xmesh 
  *****************************************************/
  for (n1=0; n1<loop_Ngrid1; n1++){

    N_loop_Ngrid1 = n1;

    if (flag_save_memory_3d==1){
      for (i=1; i<=3; i++){
        Grid_Origin[i] = Grid_Origin_N[n1][i];
      }
    }
  
    /*****************************************************
                     Set_Grid
    *****************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("start Set_Grid\n");
    }
#endif

    Set_Grid();

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("end Set_Grid\n");
      printf("\n");fflush(stdout);
    }
#endif

    MPI_Barrier(MPI_COMM_WORLD);
    /*****************************************************
                   set ddOrbitals_Grid
    *****************************************************/
#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("start Set_ddOrbitals_Grid\n");
    }
#endif

    Set_ddOrbitals_Grid();

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("end Set_ddOrbitals_Grid\n");
      printf("\n");fflush(stdout);
    }
#endif

    /*****************************************************
                    set Density Grid
    *****************************************************/

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("start Set_Density_Grid\n");
    }
#endif

    //Set_Density_Grid(DM[0]);
    Set_Density_Grid_fuku(DM[0]);
    Set_dDensity_Grid(DM[0]);
    Set_ddDensity_Grid(DM[0]);
    Set_dd2Density_Grid(DM[0]);

    if ((Solver==4)||(SpinP_switch>1)){
      Set_Density_Grid_fuku_i(iDM[0]);
      Set_dDensity_Grid_i(iDM[0]);
      Set_ddDensity_Grid_i(iDM[0]);
      Set_dd2Density_Grid_i(iDM[0]);
    }

#ifdef DEBUG
    MPI_Barrier(MPI_COMM_WORLD);
    if (myid==Host_ID){  
      printf("end Set_Density_Grid\n");
      printf("\n");fflush(stdout);
    }
#endif

    /*****************************************************
             Calculate local physical quantities 
    *****************************************************/
    sprintf(DirnameLPQ,"LPQ");

    if (myid==Host_ID){  
      printf("out_LPQ step %d / %d,    ", N_loop_Ngrid1+1, loop_Ngrid1);
      printf("Grid_Origin: %f\n", Grid_Origin[1]); fflush(stdout);
    }
    out_LPQ();

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
    if (myid==Host_ID){  
      printf("start Free_Arrays_in_Set_Grid\n");fflush(stdout);
    }
#endif

    Free_Arrays_in_Set_Grid();

#ifdef DEBUG
    if (myid==Host_ID){  
      printf("end Free_Arrays_in_Set_Grid\n");fflush(stdout);
      printf("\n");fflush(stdout);
    }
#endif

  } //loop for xmesh      
  
#ifdef DEBUG
  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){  
    printf("start Free_Arrays\n");fflush(stdout);
    printf("\n");fflush(stdout);
  }
#endif

  Free_Arrays();

#ifdef DEBUG
  if (myid==Host_ID){  
    printf("end Free_Arrays\n");fflush(stdout);
    printf("\n");fflush(stdout);
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  if (myid==Host_ID){
    printf("\nThe calculation was normally finished.\n");fflush(stdout);
  }

  MPI_Finalize();

  return 0;
}
