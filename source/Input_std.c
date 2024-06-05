#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "Inputtools.h"
#include "mpi.h"
#include <omp.h>
#include "flpq.h"



#define MAXBUF 1024


void Input_std(char *file)
{
  FILE *fp,*fp_check;
  int i,j,k,itmp;
  int po=0;  /* error count */
  double r_vec[40];
  int i_vec[40],i_vec2[40];
  char *s_vec[40];
  double tmp;
  char buf[MAXBUF];
  char file_check[YOUSO10];
  int numprocs,myid;
  int numprocs1;
  int grid_vectors_unit;
  int dens_origin_unit;

  MPI_Comm_size(MPI_COMM_WORLD1,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD1,&myid);

  if (myid==Host_ID){  
    printf("*******************************************************\n"); 
    printf("       read the input file and initializing            \n");
    printf("*******************************************************\n\n"); 
  }

  /****************************************************
                       open a file
  ****************************************************/

  if (input_open(file)==0){
    MPI_Finalize(); 
    exit(0);
  }

  Num_DM_orb=0;
  Num_DM_spin=0;
  Num_DM_k=0;

  input_int("Num_DM_orb",&Num_DM_orb,0);
  input_int("Num_DM_spin",&Num_DM_spin,-1);
  input_int("Num_DM_k",&Num_DM_k,0);
  if(Num_DM_orb>0){
    flag_DMmu=1;
    if(Num_DM_spin==-1){
      printf("You have to specify Num_DM_spin.\n");
      exit(0);
    }
  }
  if(Num_DM_k>0){
    flag_DMk=1;
  }

  input_logical("DMmu_range", &flag_DMmu_range,0);

  //input_logical("DM_orb", &flag_DMmu,0);
  //if(flag_DMmu==1){
  //  input_int("Num_DM_orb",&Num_DM_orb,0);
  //}
  //input_logical("DM_k", &flag_DMk,0);
  //if(flag_DMk==1){
  //  input_int("Num_DM_k",&Num_DM_k,0);
  //}
  //input_logical("DM_k_orb", &flag_DMkmu,0);
  //if(flag_DMkmu==1){
  //  input_int("Num_DM_k",&Num_DM_k,0);
  //  input_int("Num_DM_orb",&Num_DM_orb,0);
  //}

  input_logical("save_memory_3d", &flag_save_memory_3d,0);

  s_vec[0]="Ang"; s_vec[1]="AU";
  i_vec[0]=0;  i_vec[1]=1;
  input_string2int("dens_origin.unit",&dens_origin_unit,2,s_vec,i_vec);

  r_vec[0]=1.0e+9; r_vec[1]=1.0e+9; r_vec[2]=1.0e+9;
  input_doublev("dens_origin",3,scf_fixed_dens_origin,r_vec);

  /* Ang to AU */
  if (dens_origin_unit==0){
    scf_fixed_dens_origin[0] = scf_fixed_dens_origin[0]/BohrR;
    scf_fixed_dens_origin[1] = scf_fixed_dens_origin[1]/BohrR;
    scf_fixed_dens_origin[2] = scf_fixed_dens_origin[2]/BohrR;
  }

  i_vec[0]=0; i_vec[1]=0, i_vec[2]=0;
  input_intv("Dens_Ngrid",3,Dens_Ngrid,i_vec);

  s_vec[0]="Ang"; s_vec[1]="AU";
  i_vec[0]=0;  i_vec[1]=1;
  input_string2int("Grid_vectors.Unit",&grid_vectors_unit,2,s_vec,i_vec);

  if (fp=input_find("<Grid_vectors")) {
    flag_grid_vec = 1;

    for (i=1; i<=3; i++){
      fscanf(fp,"%lf %lf %lf",&grid_vec[i][1],&grid_vec[i][2],&grid_vec[i][3]);
    }
    if ( ! input_last("Grid_vectors>") ) {
      /* format error */
      printf("Format error for Grid_vectors\n");
      po++;
    }

    /* Ang to AU */
    if (grid_vectors_unit==0){
      for (i=1; i<=3; i++){
        grid_vec[i][1] = grid_vec[i][1]/BohrR;
        grid_vec[i][2] = grid_vec[i][2]/BohrR;
        grid_vec[i][3] = grid_vec[i][3]/BohrR;
      }
    }
  }
  else{
    flag_grid_vec = 0;
  }
 

  /****************************************************
                       input_close
  ****************************************************/

  input_close();

  if (po>0 || input_errorCount()>0) {
    printf("errors in the inputfile\n");
    MPI_Finalize();
    exit(1);
  } 

  /****************************************************
                   print out to std
  ****************************************************/

  if (myid==Host_ID){  
    printf("\n\n<Input_std>  Your input file was normally read.\n");
  }

}



