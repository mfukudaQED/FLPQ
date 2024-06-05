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
#include <time.h>
#include "flpq.h"
#include "mpi.h"
#include <unistd.h>
#include <omp.h>
#include "Local_physical_quantities.h"
#include <sys/stat.h>

static void calc_LPQ(char fileLPQ[YOUSO10], int norder, int numLPQ, int myid);
static void Print_LPQ_3D(FILE *fp, char fext[], int norder, int numLPQ);
static void Print_CubeTitle(FILE *fp);
static void Print_CubeData(FILE *fp, char fext[], int numLPQ);
static void Merge_Files(FILE *fp, char fext[], int numLPQ);
static void Print_Density_Grid();

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
     calc_LPQ(".energy",                     0, 7, myid);   /* added by takamatsu 19.February.2020 */
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



void calc_LPQ(char fileLPQ[YOUSO10], int norder, int numLPQ, int myid)
{
  char fname[YOUSO10];
  char fileLPQext[YOUSO10];
  FILE *fp;
  int fd;

    /* for QEDynamics, LPQ.*.dat */
    sprintf(fileLPQext,"%s.dat",fileLPQ);
    sprintf(fname,"%s/LPQ%s_%i",DirnameLPQ,fileLPQext,myid);
    if ((fp = fopen(fname,"w")) != NULL){
      Print_LPQ_3D(fp, fileLPQext, norder, numLPQ);
      fflush(fp);
      fd = fileno(fp); 
      fsync(fd);
      fclose(fp);

      MPI_Barrier(mpi_comm_level1);
      if (myid==Host_ID){
        Merge_Files(fp, fileLPQext, numLPQ);
      }
    }
    else{
      printf("Failure of saving the %s\n",fname);
    }

    if (flag_save_memory_3d==0){
      if(norder==0){ /* scalar */
        if((Ngrid1!=1)&&(Ngrid2!=1)&&(Ngrid3!=1)){
          /* for QEDynamics, LPQ.*.cube */
          sprintf(fileLPQext,"%s.cube",fileLPQ);
          sprintf(fname,"%s/LPQ%s_%i",DirnameLPQ,fileLPQext,myid);
          if ((fp = fopen(fname,"w")) != NULL){
            if (myid==Host_ID){
              Print_CubeTitle(fp);
            }

            Print_CubeData(fp, fileLPQext, numLPQ);

            fflush(fp);
            fd = fileno(fp); 
            fsync(fd);
            fclose(fp);

            if (myid==Host_ID){
              Merge_Files(fp, fileLPQext, numLPQ);
            }
          }
          else{
            printf("Failure of saving the %s\n",fname);
          }
        }
      }
    }

    if ((myid==Host_ID)&&(N_loop_Ngrid1==0)){
      sprintf(fname,"%s/header.cube",DirnameLPQ);
      if ((fp = fopen(fname,"w")) != NULL){
        Print_CubeTitle(fp);
        fflush(fp);
        fclose(fp);
      }

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

static void Merge_Files(FILE *fp, char fext[], int numLPQ)
{
  int numprocs,ID,digit,c;
  char operate[500],operate1[500],operate2[500];
  FILE *fp1,*fp2;

  /****************************************************
                   merge files
  ****************************************************/

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  digit = (int)log10(numprocs) + 1;


  //sprintf(operate,"cat %s/LPQ%s_0 > %s/LPQ%s",
  //        DirnameLPQ, fext, DirnameLPQ,fext);
  //system(operate);

  //for (ID=1; ID<numprocs; ID++){
  //  sprintf(operate,"cat %s/LPQ%s_%i >> %s/LPQ%s",
  //          DirnameLPQ,fext,ID, DirnameLPQ,fext);
  //  system(operate);
  //}

  //for (ID=0; ID<numprocs; ID++){
  //  sprintf(operate,"rm %s/LPQ%s_%i",DirnameLPQ, fext,ID);
  //  system(operate);
  //}


  sprintf(operate1,"%s/LPQ%s",DirnameLPQ,fext);
  //fp1 = fopen(operate1, "r");   

  //if (fp1!=NULL){
  //  fclose(fp1); 
  //  remove(operate1);
  //}

  /* merge all the fraction files */

  for (ID=0; ID<numprocs; ID++){

    sprintf(operate1,"%s/LPQ%s",DirnameLPQ,fext);
    fp1 = fopen(operate1, "a");   
    fseek(fp1,0,SEEK_END);

    sprintf(operate2,"%s/LPQ%s_%i",DirnameLPQ,fext,ID);
    fp2 = fopen(operate2,"r");

    if (fp2!=NULL){
      for (c=getc(fp2); c!=EOF; c=getc(fp2))  putc(c,fp1); 
      fclose(fp2); 
    }

    fclose(fp1); 
  
  }

  for (ID=0; ID<numprocs; ID++){
    sprintf(operate,"%s/LPQ%s_%i",DirnameLPQ,fext,ID);
    remove(operate);
  }

}
#endif



