/**********************************************************************
  out_LPQ.c:

     out_LPQ.c is a subrutine to output values of local physical quantities.

  Log of out_LPQ.c:

     01/Jul/2016  Released by M.Fukuda
     18/Jun/2019  Modified by M.Fukuda
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <ctype.h>
#include <time.h>
#include "flpq.h"
#include "mpi.h"
#include <unistd.h>
#include <omp.h>
#include "Local_physical_quantities.h"
#include "FLPQ_Quantity_Selector.h"
#include "FLPQ_HDF5.h"
#include <sys/stat.h>
#include <errno.h>

typedef struct {
  char LPQNAME[50];
  int dimension;
  int index;
  int SpinP_switch;
} Struct_LPQ;

static void calc_LPQ(const char *fileLPQ, int norder, int numLPQ, int myid);
static void Print_LPQ_3D(FILE *fp, char fext[], int norder, int numLPQ);
static void Print_CubeTitle(FILE *fp);
static void Print_CubeData(FILE *fp, char fext[], int numLPQ);
static void Merge_Files(const char *fext);
static void Print_Density_Grid();
static void flpq_quantity_error(const char *msg, const char *quantity)
{
  int me;
  MPI_Comm_rank(mpi_comm_level1,&me);
  fprintf(stderr,"[FLPQ rank %d] Output.Quantities: %s %s\n",
          me,msg,quantity?quantity:"");
  fflush(stderr);
  MPI_Abort(mpi_comm_level1,EXIT_FAILURE);
}


/* added by fukuda 11.May.2015 */
void out_LPQ()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fileLPQ[YOUSO10];
  int numprocs,myid,ID;
  char operate[100];
  int dircheck;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if ((myid==Host_ID)&&(N_loop_Ngrid1==0)){
    //sprintf(operate,"mkdir %s",DirnameLPQ);
    mkdir(DirnameLPQ,0755);
  }
  MPI_Barrier(mpi_comm_level1);

  /****************************************************
               Local physical quantities 
  ****************************************************/
     calc_LPQ(".electron_density",           0, 0, myid);
     calc_LPQ(".charge_density",             0, 1, myid);
     if(SpinP_switch==3)
     calc_LPQ(".zeta_potential",             0, 2, myid);
     calc_LPQ(".kinetic_energy_density",     0, 3, myid);
//     //calc_LPQ(".energy_density",         0, 4, myid);
     calc_LPQ(".charge_density_alpha",       0, 5, myid);   /* added by yoshida 23.December.2019 */
     calc_LPQ(".charge_density_beta",        0, 6, myid);   /* added by yoshida 23.December.2019 */
     calc_LPQ(".energy_density_tau",         0, 7, myid);   /* added by takamatsu 19.February.2020 */
     calc_LPQ(".chemical_potential",         0, 8, myid);   /* added by takamatsu 19.February.2020 */     
     calc_LPQ(".charge_current",             1, 0, myid);
     calc_LPQ(".kinetic_momentum",           1, 1, myid);
     calc_LPQ(".spin_vorticity",             1, 2, myid);
     calc_LPQ(".spin_angular_momentum",      1, 3, myid);
//     //calc_LPQ(".spin_torque",            1, 4, myid);
     if(SpinP_switch==3)
     calc_LPQ(".zeta_force",                 1, 5, myid);
     if(SpinP_switch==3)
     calc_LPQ(".pau_dd",                     1, 6, myid);
//     //calc_LPQ(".tension",                1, 7, myid);
//     //calc_LPQ(".lorentz_force",          1, 8, myid);
     calc_LPQ(".gradient_electron_density",  1, 9, myid);

     if(SpinP_switch==3){
       calc_LPQ(".spin_z_current",             1, 10, myid);  /* added by kuroda 7.June.2019 */
       calc_LPQ(".charge_current_alpha",       1, 11, myid);  /* added by kuroda 5.December.2019 */
       calc_LPQ(".charge_current_beta",        1, 12, myid);  /* added by kuroda 5.December.2019 */
       calc_LPQ(".velocity",                   1, 13, myid);  /* added by kuroda 31.January.2020 */
       calc_LPQ(".torq_zeta",                  1, 14, myid);  /* added by shimada 27.December.2019 */ 
       calc_LPQ(".soc_part1",                  1, 15, myid);  /*added by shimada September.2020*/
       calc_LPQ(".kinetic_momentum_alpha",     1, 16, myid);
       calc_LPQ(".kinetic_momentum_beta",      1, 17, myid);
       calc_LPQ(".soc_part2",                  1, 18, myid);  /*added by shimada September.2020*/
    }
    
     calc_LPQ(".stress",                     2, 0, myid);
     calc_LPQ(".stress_S",                   2, 1, myid);
     if(SpinP_switch==3)
     calc_LPQ(".stress_A",                   2, 2, myid);
     calc_LPQ(".stress_diag",                2, 3, myid);
     calc_LPQ(".dd_electron_density",        2, 4, myid);
     calc_LPQ(".stress_nonrel",              2, 5, myid);   /* added by yoshida 7.January.2020 */
     calc_LPQ(".stress_nonrel_alpha",        2, 6, myid);   /* added by yoshida 7.January.2020 */ 
     calc_LPQ(".stress_nonrel_beta",         2, 7, myid);   /* added by yoshida 7.January.2020 */
}

/* modified by fukuda 08.July.2024 */
void out_LPQ2()
{
  int ct_AN,spe,i1,i2,i3,GN,i,j,k;
  char fileLPQ[YOUSO10];
  int numprocs,myid,ID;
  char operate[100];
  int dircheck;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if ((myid==Host_ID)&&(N_loop_Ngrid1==0)){
    //sprintf(operate,"mkdir %s",DirnameLPQ);
    mkdir(DirnameLPQ,0755);
  }
  MPI_Barrier(mpi_comm_level1);

  /****************************************************
               Local physical quantities 
  ****************************************************/
  Struct_LPQ dataLPQ[] = {
    {".electron_density",           0, 0, 0},
    //{".charge_density",             0, 1, 0},
    {".zeta_potential",             0, 2, 3},
    {".kinetic_energy_density",     0, 3, 0},
    //{".charge_density_alpha",       0, 5, 0},
    //{".charge_density_beta",        0, 6, 0},
    {".energy_density_tau",         0, 7, 0},
    {".chemical_potential",         0, 8, 0},

    {".charge_current",             1, 0, 0},
    {".kinetic_momentum",           1, 1, 0},
    {".spin_vorticity",             1, 2, 0},
    {".spin_angular_momentum",      1, 3, 0},
    {".zeta_force",                 1, 5, 3},
    //{".pau_dd",                     1, 6, 3},
    //{".gradient_electron_density",  1, 9, 0},

    //{".spin_z_current",             1, 10, 3},
    //{".charge_current_alpha",       1, 11, 3},
    //{".charge_current_beta",        1, 12, 3},
    //{".velocity",                   1, 13, 3},
    //{".torq_zeta",                  1, 14, 3},
    //{".soc_part1",                  1, 15, 3},
    //{".kinetic_momentum_alpha",     1, 16, 3},
    //{".kinetic_momentum_beta",      1, 17, 3},
    //{".soc_part2",                  1, 18, 3},

    //{".stress",                     2, 0, 0},
    //{".stress_S",                   2, 1, 0},
    //{".stress_A",                   2, 2, 3},
    {".stress_diag",                2, 3, 0},
    //{".dd_electron_density",        2, 4, 0},
    //{".stress_nonrel",              2, 5, 0},
    //{".stress_nonrel_alpha",        2, 6, 0},
    //{".stress_nonrel_beta",         2, 7, 0}
  };

  /* Ordinary defaults remain identical to F3; decomposition remains opt-in. */
  const Struct_LPQ optional_stress[] = {
    {".stress",2,0,0}, {".stress_S",2,1,0}, {".stress_A",2,2,0}
  };
  int size_dataLPQ = (int)(sizeof(dataLPQ) / sizeof(dataLPQ[0]));
  int selected[32] = {0}, ntypes=size_dataLPQ+3;
  const Struct_LPQ *item[32];
  FLPQ_Quantity_Descriptor descriptors[32];
  char selection_error[256];
  if(ntypes>32)flpq_quantity_error("too many registered quantities",NULL);
  for(i=0;i<size_dataLPQ;i++) item[i]=&dataLPQ[i];
  for(i=0;i<3;i++) item[size_dataLPQ+i]=&optional_stress[i];
  for(i=0;i<ntypes;i++){
    descriptors[i].name=item[i]->LPQNAME+1;
    descriptors[i].dimension=item[i]->dimension;
    descriptors[i].noncollinear_only=item[i]->SpinP_switch==3;
    descriptors[i].optional_stress=i>=size_dataLPQ;
  }
  if(flpq_select_output_quantities(FLPQ_Output_Quantities,descriptors,ntypes,
                                   SpinP_switch,FLPQ_Output_Stress_Decomposition,
                                   selected,selection_error,sizeof(selection_error))!=0)
    flpq_quantity_error(selection_error,NULL);
  /* Fail before generating any selected HDF5 file, not partway through.
   * The .partial file is deliberately never deleted automatically.
   * Explicit replace only swaps the final file after successful HDF5 close.
   */
  if (myid==Host_ID && N_loop_Ngrid1==0 && FLPQ_Output_Format!=0){
    for(i=0;i<ntypes;i++) if(selected[i]){
      char final[YOUSO10],partial[YOUSO10];
      struct stat st;
      if(snprintf(final,sizeof(final),"%s/LPQ%s.h5",DirnameLPQ,
                  item[i]->LPQNAME)>=(int)sizeof(final) ||
         snprintf(partial,sizeof(partial),"%s/LPQ%s.partial.h5",
                  DirnameLPQ,item[i]->LPQNAME)>=(int)sizeof(partial))
        flpq_quantity_error("output path too long",item[i]->LPQNAME);
      if(lstat(partial,&st)==0)
        flpq_quantity_error("partial HDF5 exists; intervention needed",partial);
      if(lstat(final,&st)==0){
        if(FLPQ_HDF5_OnExisting==0)
          flpq_quantity_error("HDF5 exists, refusing overwrite",final);
        if(!S_ISREG(st.st_mode))
          flpq_quantity_error("replacement target is not a regular file",final);
      }
    }
  }
  MPI_Barrier(mpi_comm_level1);
  for(i=0;i<ntypes;i++) if(selected[i])
    calc_LPQ(item[i]->LPQNAME,item[i]->dimension,item[i]->index,myid);

}

static void flpq_output_fail(const char *action, const char *path)
{
  int rank;
  MPI_Comm_rank(mpi_comm_level1, &rank);
  fprintf(stderr, "[FLPQ rank %d] %s %s: %s\n",
          rank, action, path, strerror(errno));
  fflush(stderr);
  MPI_Abort(mpi_comm_level1, EXIT_FAILURE);
}

static void flpq_output_close(FILE *fp, const char *path)
{
  if (fclose(fp) != 0) flpq_output_fail("cannot close", path);
}

static void calc_LPQ(const char *fileLPQ, int norder, int numLPQ, int myid)
{
  char fname[YOUSO10], fileLPQext[YOUSO10];
  FILE *fp;

  if (FLPQ_Output_Format != 0)
    FLPQ_HDF5_write(fileLPQ, norder, numLPQ);
  if (FLPQ_Output_Format == 1) return;

  if (snprintf(fileLPQext, sizeof(fileLPQext), "%s.dat", fileLPQ)
      >= (int)sizeof(fileLPQext) ||
      snprintf(fname, sizeof(fname), "%s/LPQ%s_%d",
               DirnameLPQ, fileLPQext, myid) >= (int)sizeof(fname)) {
    errno = ENAMETOOLONG;
    flpq_output_fail("output path too long", fileLPQ);
  }
  fp = fopen(fname, "w");
  if (fp == NULL) flpq_output_fail("cannot open", fname);
  Print_LPQ_3D(fp, fileLPQext, norder, numLPQ);
  flpq_output_close(fp, fname);

  /* Every fragment must be closed before rank zero starts merging. */
  MPI_Barrier(mpi_comm_level1);
  if (myid == Host_ID) Merge_Files(fileLPQext);
  MPI_Barrier(mpi_comm_level1);

  if (flag_save_memory_3d == 0 && norder == 0 &&
      Ngrid1 != 1 && Ngrid2 != 1 && Ngrid3 != 1) {
    if (snprintf(fileLPQext, sizeof(fileLPQext), "%s.cube", fileLPQ)
        >= (int)sizeof(fileLPQext) ||
        snprintf(fname, sizeof(fname), "%s/LPQ%s_%d",
                 DirnameLPQ, fileLPQext, myid) >= (int)sizeof(fname)) {
      errno = ENAMETOOLONG;
      flpq_output_fail("output path too long", fileLPQ);
    }
    fp = fopen(fname, "w");
    if (fp == NULL) flpq_output_fail("cannot open", fname);
    if (myid == Host_ID) Print_CubeTitle(fp);
    Print_CubeData(fp, fileLPQext, numLPQ);
    flpq_output_close(fp, fname);

    MPI_Barrier(mpi_comm_level1);
    if (myid == Host_ID) Merge_Files(fileLPQext);
    MPI_Barrier(mpi_comm_level1);
  }

  if (myid == Host_ID && N_loop_Ngrid1 == 0) {
    if (snprintf(fname, sizeof(fname), "%s/header.cube", DirnameLPQ)
        >= (int)sizeof(fname)) {
      errno = ENAMETOOLONG;
      flpq_output_fail("output path too long", DirnameLPQ);
    }
    fp = fopen(fname, "w");
    if (fp == NULL) flpq_output_fail("cannot open", fname);
    Print_CubeTitle(fp);
    flpq_output_close(fp, fname);
  }
}

/* added by fukuda 21.May.2023 */
static void Print_Density_Grid()
{
  FILE *fp;
  int GN,spin;
  //double x,y,z;
  double Cxyz[4];
  int numprocs,myid,ID;
  char operate[300];
  char fname[500];
  int SpinP_Grid_size;
  int i;


  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                    output data 
  ****************************************************/
  sprintf(fname,"Density_Grid_%i.bin",myid);
  fp = fopen(fname,"wb");

  if (myid==Host_ID){
    fwrite(&SpinP_switch, sizeof(int),1,fp);
    fwrite(&SpinP_Grid_size, sizeof(int),1,fp);
    fwrite(&Ngrid1, sizeof(int),1,fp);
    fwrite(&Ngrid2, sizeof(int),1,fp);
    fwrite(&Ngrid3, sizeof(int),1,fp);
    fwrite(&Grid_Origin[1], sizeof(double),3,fp);
    fwrite(&gtv[1][1], sizeof(double),3,fp);
    fwrite(&gtv[2][1], sizeof(double),3,fp);
    fwrite(&gtv[3][1], sizeof(double),3,fp);
    //for (GN=0; GN<Ngrid1*Ngrid2*Ngrid3; GN++){
    //  Get_Grid_XYZ(GN,Cxyz);
    //  //x = Cxyz[1];
    //  //y = Cxyz[2];
    //  //z = Cxyz[3];
    //  fwrite(&Cxyz[1], sizeof(double),3,fp);
    //}
  }

  if (SpinP_switch==3){
       SpinP_Grid_size = SpinP_switch + 2;
  }
  else{
       SpinP_Grid_size = SpinP_switch;
  }

  for (i=0; i<My_NumGridB_AB; i++){
 
    fwrite(&Density_Grid_B[i], sizeof(double),(SpinP_switch+1),fp);
    fwrite(&Density_Grid_B_i[i], sizeof(double),2,fp);

    fwrite(&dDensity_Grid_B[i], sizeof(double),3*(SpinP_Grid_size+1),fp);
    fwrite(&dDensity_Grid_B_i[i], sizeof(double),3*2,fp);

    fwrite(&ddDensity_Grid_B[i], sizeof(double),7*(SpinP_Grid_size+1),fp);
    fwrite(&ddDensity_Grid_B_i[i], sizeof(double),7*2,fp);

    fwrite(&dd2Density_Grid_B[i], sizeof(double),7*(SpinP_Grid_size+1),fp);
    fwrite(&dd2Density_Grid_B_i[i], sizeof(double),7*2,fp);

  }
 
  //fwrite(&Density_Grid_B, sizeof(double),My_NumGridB_AB*(SpinP_switch+1),fp);
  //fwrite(&Density_Grid_B_i, sizeof(double),My_NumGridB_AB*2,fp);

  //fwrite(&dDensity_Grid_B, sizeof(double),My_NumGridB_AB*3*(SpinP_Grid_size+1),fp);
  //fwrite(&dDensity_Grid_B_i, sizeof(double),My_NumGridB_AB*3*2,fp);

  //fwrite(&ddDensity_Grid_B, sizeof(double),My_NumGridB_AB*7*(SpinP_Grid_size+1),fp);
  //fwrite(&ddDensity_Grid_B_i, sizeof(double),My_NumGridB_AB*7*2,fp);

  //fwrite(&dd2Density_Grid_B, sizeof(double),My_NumGridB_AB*7*(SpinP_Grid_size+1),fp);
  //fwrite(&dd2Density_Grid_B_i, sizeof(double),My_NumGridB_AB*7*2,fp);

  /****************************************************
  fclose(fp);
  ****************************************************/

  fclose(fp);

  /****************************************************
                   merge files
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);

  if (myid==Host_ID){

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"cat %s/Density_Grid_%i.bin >> %s/Density_Grid.bin",
              DirnameLPQ,ID,DirnameLPQ);
      system(operate);
    }

    for (ID=0; ID<numprocs; ID++){
      sprintf(operate,"rm %s/Density_Grid_%i.bin",DirnameLPQ,ID);
      system(operate);
    }
  }
}

/* added by fukuda 2.May.2015 */
static void Print_LPQ_3D(FILE *fp, char fext[], int norder, int numLPQ)
{
  int i,j,k,i1,i2,i3,c,GridNum;
  int GN_AB,n1,n2,n3,interval;
  int BN_AB,N2D,GNs,GN;
  double x,y,z,vx;
  double Cxyz[4];
  int numprocs,myid,ID;
  char operate[300];


  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                    output data 
  ****************************************************/

  /* for QEDynamics */ 

  if ((myid==Host_ID)&&(N_loop_Ngrid1==0)){
    if (flag_save_memory_3d==0){
      fprintf(fp,"# Ngrid %i %i %i\n",Ngrid1,Ngrid2,Ngrid3);
    }
    else {
      fprintf(fp,"# Ngrid %i %i %i\n",loop_Ngrid1,Ngrid2,Ngrid3);
    }
    if (norder==0) fprintf(fp,"# x, y, z, scalar (unit : [bohr], [a.u.])\n");
    else if (norder==1) fprintf(fp,"# x, y, z, vector (unit : [bohr], [a.u.])\n");
    else if (norder==2) fprintf(fp,"# x, y, z, tensor (unit : [bohr], [a.u.]) (xx,yx,zx,xy,yy,zy,xz,yz,zz)\n");
  }

  /****************************************************
                 fprintf scalar data
  ****************************************************/

  N2D = Ngrid1*Ngrid2;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*Ngrid3;

  for (BN_AB=0; BN_AB<My_NumGridB_AB; BN_AB++){
   
    GN_AB = BN_AB + GNs;
    n1 = GN_AB/(Ngrid2*Ngrid3);
    n2 = (GN_AB - n1*(Ngrid2*Ngrid3))/Ngrid3;
    n3 = GN_AB - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;


    GN = n1*Ngrid2*Ngrid3 + n2*Ngrid3 + n3; 

    Get_Grid_XYZ(GN,Cxyz);
    x = Cxyz[1];
    y = Cxyz[2];
    z = Cxyz[3];

    //fprintf(fp,"%18.8E %18.8E %18.8E", x,y,z);
    fprintf(fp,"%14.4E %14.4E %14.4E", x,y,z);

    /* calculate values of local physical quantity */
    if (norder==0) scalar_density(numLPQ, fp, BN_AB, 0);
    else if (norder==1) vector_density(numLPQ, fp, BN_AB);
    else if (norder==2) tensor_density(numLPQ, fp, BN_AB);

    if (Ngrid3==1) {
      if (GN%(Ngrid2*Ngrid3)==(Ngrid2*Ngrid3-1)) {
        fprintf(fp,"\n");
      }
    }
    else if (GN%Ngrid3==(Ngrid3-1)) {
       fprintf(fp,"\n");
    }

  }

  ///****************************************************
  //                 merge files
  //****************************************************/

  //MPI_Barrier(mpi_comm_level1);

  //if (myid==Host_ID){

  //  sprintf(operate,"cat %s/LPQ%s0 > %s/LPQ%s",
  //          DirnameLPQ, fext, DirnameLPQ,fext);
  //  system(operate);

  //  for (ID=1; ID<numprocs; ID++){
  //    sprintf(operate,"cat %s/LPQ%s%i >> %s/LPQ%s",
  //            DirnameLPQ,fext,ID, DirnameLPQ,fext);
  //    system(operate);
  //  }

  //  for (ID=0; ID<numprocs; ID++){
  //    sprintf(operate,"rm %s/LPQ%s%i",DirnameLPQ, fext,ID);
  //    system(operate);
  //  }
  //}
}




#if 1
/* added by fukuda 17.June.2019 */
static void Print_CubeTitle(FILE *fp)
{
  int ct_AN;
  int spe; 

  fprintf(fp," SYS1\n SYS1\n");

  if (flag_save_memory_3d==0){
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        atomnum,Grid_Origin[1],Grid_Origin[2],Grid_Origin[3]);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        Ngrid1,gtv[1][1],gtv[1][2],gtv[1][3]);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        Ngrid2,gtv[2][1],gtv[2][2],gtv[2][3]);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        Ngrid3,gtv[3][1],gtv[3][2],gtv[3][3]);
  }
  else{
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        atomnum,Grid_Origin_N[0][1],Grid_Origin_N[0][2],Grid_Origin_N[0][3]);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        loop_Ngrid1,tv[1][1]/(double)loop_Ngrid1,tv[1][2]/(double)loop_Ngrid1,tv[1][3]/(double)loop_Ngrid1);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        Ngrid2,tv[2][1]/(double)Ngrid2,tv[2][2]/(double)Ngrid2,tv[2][3]/(double)Ngrid2);
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf\n",
        Ngrid3,tv[3][1]/(double)Ngrid3,tv[3][2]/(double)Ngrid3,tv[3][3]/(double)Ngrid3);
  }

  for (ct_AN=1; ct_AN<=atomnum; ct_AN++){
    spe = WhatSpecies[ct_AN];
    fprintf(fp,"%5d%12.6lf%12.6lf%12.6lf%12.6lf\n",
	    Spe_WhatAtom[spe],
	    Spe_Core_Charge[spe]-InitN_USpin[ct_AN]-InitN_DSpin[ct_AN],
	    Gxyz[ct_AN][1],Gxyz[ct_AN][2],Gxyz[ct_AN][3]);
  }

}


/* added by fukuda 17.June.2019 */
static void Print_CubeData(FILE *fp, char fext[], int numLPQ)
{
  int i,j,k,i1,i2,i3,c,GridNum;
  int GN_AB,n1,n2,n3,interval;
  int BN_AB,N2D,GNs,GN;
  double x,y,z,vx;
  double Cxyz[4];
  int numprocs,myid,ID;
  char operate[300];
  int fd;


  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
                 fprintf scalar data
  ****************************************************/

  N2D = Ngrid1*Ngrid2;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*Ngrid3;

  //printf("myid=%d, %d\n",myid,My_NumGridB_AB);fflush(stdout);

  for (BN_AB=0; BN_AB<My_NumGridB_AB; BN_AB+=Ngrid3){
    for (n3=0; n3<Ngrid3; n3++){
   
      //GN_AB = BN_AB + GNs;
      //n1 = GN_AB/(Ngrid2*Ngrid3);
      //n2 = (GN_AB - n1*(Ngrid2*Ngrid3))/Ngrid3;
      ////n3 = GN_AB - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;

      /* calculate values of local physical quantity */
      scalar_density(numLPQ, fp, BN_AB+n3, 1); /* the last number "1" means that flag_cube is ON. */

      if ((n3+1)%6==0) { fprintf(fp,"\n"); }
    }
    /* avoid double \n\n when Ngrid3%6 == 0  */
    if (Ngrid3%6!=0) fprintf(fp,"\n");

  }

  //fd = fileno(fp); 
  //fsync(fd);
  //fflush(stdout);
  //fclose(fp);

  /****************************************************
                   merge files
  ****************************************************/

  //MPI_Barrier(mpi_comm_level1);

  //if (myid==Host_ID){

  //  sprintf(operate,"cat %s/LPQ%s0 > %s/LPQ%s",
  //          DirnameLPQ, fext, DirnameLPQ,fext);
  //  system(operate);
  //  fflush(stdout);

  //  for (ID=1; ID<numprocs; ID++){
  //    sprintf(operate,"cat %s/LPQ%s%i >> %s/LPQ%s",
  //            DirnameLPQ,fext,ID, DirnameLPQ,fext);
  //    system(operate);
  //    fflush(stdout);
  //  }

  //  //for (ID=0; ID<numprocs; ID++){
  //  //  sprintf(operate,"rm %s/LPQ%s%i",DirnameLPQ, fext,ID);
  //  //  system(operate);
  //  //}
  //}
}

static void Merge_Files(const char *fext)
{
  int numprocs, rank;
  char output_path[YOUSO10], fragment_path[YOUSO10];
  char buffer[65536];
  FILE *output, *fragment;
  size_t count;

  MPI_Comm_size(mpi_comm_level1, &numprocs);
  if (snprintf(output_path, sizeof(output_path), "%s/LPQ%s",
               DirnameLPQ, fext) >= (int)sizeof(output_path)) {
    errno = ENAMETOOLONG;
    flpq_output_fail("output path too long", fext);
  }

  /* The first slice replaces old output; later slices are appended. */
  output = fopen(output_path, N_loop_Ngrid1 == 0 ? "wb" : "ab");
  if (output == NULL) flpq_output_fail("cannot open", output_path);

  for (rank = 0; rank < numprocs; rank++) {
    if (snprintf(fragment_path, sizeof(fragment_path), "%s/LPQ%s_%d",
                 DirnameLPQ, fext, rank) >= (int)sizeof(fragment_path)) {
      errno = ENAMETOOLONG;
      flpq_output_fail("fragment path too long", fext);
    }
    fragment = fopen(fragment_path, "rb");
    if (fragment == NULL) flpq_output_fail("cannot open", fragment_path);

    while ((count = fread(buffer, 1, sizeof(buffer), fragment)) != 0) {
      if (fwrite(buffer, 1, count, output) != count)
        flpq_output_fail("cannot write", output_path);
    }
    if (ferror(fragment)) flpq_output_fail("cannot read", fragment_path);
    flpq_output_close(fragment, fragment_path);
  }
  flpq_output_close(output, output_path);

  /* Delete fragments only after the complete destination is closed. */
  for (rank = 0; rank < numprocs; rank++) {
    if (snprintf(fragment_path, sizeof(fragment_path), "%s/LPQ%s_%d",
                 DirnameLPQ, fext, rank) >= (int)sizeof(fragment_path)) {
      errno = ENAMETOOLONG;
      flpq_output_fail("fragment path too long", fext);
    }
    if (remove(fragment_path) != 0) flpq_output_fail("cannot remove", fragment_path);
  }
}
#endif



