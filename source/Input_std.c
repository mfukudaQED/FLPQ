#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include "Inputtools.h"
#include "mpi.h"
#include <omp.h>
#include "flpq.h"



#define MAXBUF 1024


int FLPQ_Output_Stress_Decomposition = 0;
int FLPQ_Output_Format = 0; /* legacy=0, hdf5=1, both=2 */
char FLPQ_Output_Quantities[1024] = "all";
int FLPQ_HDF5_Precision = 64; /* float64=64, float32=32 */
int FLPQ_HDF5_Compression = 1; /* none=0, gzip=1 */
int FLPQ_HDF5_GzipLevel = 4;
int FLPQ_HDF5_OnExisting = 0; /* error=0, replace=1 */
int FLPQ_HDF5_ChunkY = 16;
int FLPQ_HDF5_ChunkZ = 128;
int FLPQ_HDF5_BlockPoints = 512;
int FLPQ_HDF5_LogPerformance = 0;


/* Parse the whole value because Output.Quantities may contain commas and spaces.
 * Unlike the generic single-token input helpers, this rejects duplicate,
 * overlong, empty and truncated entries rather than silently dropping data. */
static int read_output_quantities(const char *name, char *result,
                                  size_t capacity, int rank)
{
  FILE *in = fopen(name, "r");
  char line[4096];
  int seen = 0, line_no = 0, i;
  const char *settings[] = {
    "Output.Format", "Output.HDF5.Precision",
    "Output.HDF5.Compression", "Output.HDF5.GzipLevel",
    "Output.HDF5.OnExisting", "Output.HDF5.ChunkY",
    "Output.HDF5.ChunkZ", "Output.HDF5.BlockPoints",
    "Output.HDF5.LogPerformance"
  };
  int seen_settings[9] = {0};
  if (in == NULL) {
    fprintf(stderr,"[FLPQ rank %d] cannot reopen input for Output.Quantities\n",rank);
    return 0;
  }
  while (fgets(line, sizeof(line), in)) {
    char *start=line, *end, *sep;
    size_t len=strlen(line), value_len;
    line_no++;
    if (len == sizeof(line)-1 && line[len-1]!='\n') {
      fprintf(stderr,"[FLPQ rank %d] overlong input line %d\n",rank,line_no);
      fclose(in);return 0;
    }
    while (isspace((unsigned char)*start)) start++;
    if (*start=='#' || *start=='\0') continue;
    end=start;
    while (*end && !isspace((unsigned char)*end)) end++;
    /* Generic Inputtools reads the first occurrence; F4 forbids ambiguity. */
    for(i=0;i<9;i++){
      if((size_t)(end-start)==strlen(settings[i]) &&
         strncasecmp(start,settings[i],(size_t)(end-start))==0){
        if(seen_settings[i]++){
          fprintf(stderr,"[FLPQ rank %d] duplicate %s on line %d\n",
                  rank,settings[i],line_no);
          fclose(in);return 0;
        }
      }
    }
    if ((size_t)(end-start)!=strlen("Output.Quantities") ||
        strncasecmp(start,"Output.Quantities",(size_t)(end-start))) continue;
    if (seen++) {
      fprintf(stderr,"[FLPQ rank %d] duplicate Output.Quantities on line %d\n",
              rank,line_no);
      fclose(in);return 0;
    }
    start=end;
    while (isspace((unsigned char)*start)) start++;
    sep=strchr(start,'#');if(sep)*sep='\0';
    end=start+strlen(start);
    while(end>start && isspace((unsigned char)end[-1])) end--;
    *end='\0';
    value_len=(size_t)(end-start);
    if(value_len==0 || value_len>=capacity){
      fprintf(stderr,"[FLPQ rank %d] invalid/overlong Output.Quantities on line %d\n",
              rank,line_no);
      fclose(in);return 0;
    }
    memcpy(result,start,value_len+1);
  }
  if (ferror(in)) {
    fprintf(stderr,"[FLPQ rank %d] failed reading Output.Quantities\n",rank);
    fclose(in);return 0;
  }
  fclose(in);
  return 1;
}

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
    fprintf(stderr, "[FLPQ rank %d] cannot open input file %s\n", myid, file);
    fflush(stderr);
    MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
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
      fprintf(stderr, "[FLPQ rank %d] Num_DM_spin is required with Num_DM_orb\n", myid);
      MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
    }
  }
  if(Num_DM_k>0){
    flag_DMk=1;
  }

  if (Num_DM_orb<0 || Num_DM_k<0) {
    fprintf(stderr, "[FLPQ rank %d] negative density-matrix selection\n", myid);
    MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
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
  input_logical("Output.Stress.Decomposition", &FLPQ_Output_Stress_Decomposition,0);
  /* HDF5 v1 is opt-in; legacy remains the default. */
  s_vec[0]="legacy"; s_vec[1]="hdf5"; s_vec[2]="both";
  i_vec[0]=0; i_vec[1]=1; i_vec[2]=2;
  if (input_string2int("Output.Format", &FLPQ_Output_Format, 3, s_vec, i_vec) < 0) po++;
  s_vec[0]="float64"; s_vec[1]="float32";
  i_vec[0]=64; i_vec[1]=32;
  if (input_string2int("Output.HDF5.Precision", &FLPQ_HDF5_Precision, 2, s_vec, i_vec) < 0) po++;
  /* input_string2int uses entry 0 as its default; gzip is the default. */
  s_vec[0]="gzip"; s_vec[1]="none";
  i_vec[0]=1; i_vec[1]=0;
  if (input_string2int("Output.HDF5.Compression", &FLPQ_HDF5_Compression, 2, s_vec, i_vec) < 0) po++;
  if (input_int("Output.HDF5.GzipLevel", &FLPQ_HDF5_GzipLevel, 4) < 0) po++;
  s_vec[0]="error"; s_vec[1]="replace";
  i_vec[0]=0; i_vec[1]=1;
  if (input_string2int("Output.HDF5.OnExisting", &FLPQ_HDF5_OnExisting, 2, s_vec, i_vec) < 0) po++;
  if (input_int("Output.HDF5.ChunkY", &FLPQ_HDF5_ChunkY, 16) < 0) po++;
  if (input_int("Output.HDF5.ChunkZ", &FLPQ_HDF5_ChunkZ, 128) < 0) po++;
  if (input_int("Output.HDF5.BlockPoints", &FLPQ_HDF5_BlockPoints, 512) < 0) po++;
  if (input_logical("Output.HDF5.LogPerformance", &FLPQ_HDF5_LogPerformance, 0) < 0) po++;
  if (FLPQ_HDF5_ChunkY<1 || FLPQ_HDF5_ChunkY>128 ||
      FLPQ_HDF5_ChunkZ<1 || FLPQ_HDF5_ChunkZ>512 ||
      FLPQ_HDF5_BlockPoints<16 || FLPQ_HDF5_BlockPoints>8192) {
    fprintf(stderr,"[FLPQ rank %d] invalid HDF5 I/O tuning: "
            "ChunkY=1..128 ChunkZ=1..512 BlockPoints=16..8192\n",myid);
    po++;
  }

  if(FLPQ_HDF5_GzipLevel<0 || FLPQ_HDF5_GzipLevel>9){
    fprintf(stderr,"[FLPQ rank %d] Output.HDF5.GzipLevel must be 0..9\n",myid);
    po++;
  }
  if(!read_output_quantities(file,FLPQ_Output_Quantities,
                             sizeof(FLPQ_Output_Quantities),myid)) po++;



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
    fprintf(stderr, "[FLPQ rank %d] errors in input file %s\n", myid, file);
    MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
  } 

  /****************************************************
                   print out to std
  ****************************************************/

  if (myid==Host_ID){  
    printf("\n\n<Input_std>  Your input file was normally read.\n");
  }

}



