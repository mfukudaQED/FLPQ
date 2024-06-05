/**********************************************************************
  Set_Allocate_Atom2CPU.c:

    Set_Allocate_Atom2CPU.c is a subroutine to allocate atoms to processors
    for the MPI parallel computation.

  Log of Set_Allocate_Atom2CPU.c:

     22/Nov/2001  Released by T.Ozaki
     19/May/2019  Modified by M.FUKUDA

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "flpq.h"
#include "mpi.h"
#include "lapack_prototype.h"

static void Allocation_Atoms_3D(int NL_switch);
static void Eigen_lapack(double **a, double *ko, int n, int EVmax);
static void Eigen_lapack_x(double **a, double *ko, int n0, int EVmax);
static void Eigen_lapack_d(double **a, double *ko, int n0, int EVmax);
static void Eigen_HH(double **ac, double *ko, int n, int EVmax);



int Set_Allocate_Atom2CPU(int weight_flag)
{
  double time0;
  time_t TStime,TEtime;

  time(&TStime);

  Allocation_Atoms_3D(weight_flag);
  
  time(&TEtime);
  time0 = TEtime - TStime;
  return time0;
}



#pragma optimization_level 1
void Allocation_Atoms_3D(int weight_flag)
{
  /***************************************
        allocate atoms to processes
      by Modified Recursive Bisection
  ***************************************/
  int m,i,j,k,k0,Na,np,ID,n0; 
  int myid,numprocs,numprocs0;
  int max_depth,n,depth,child;
  int **Num_Procs_in_Child;
  int **Num_Atoms_in_Child;
  int ***List_AN_in_Child;
  int *MatomN;
  int alloc_first_10;
  double t,ax,ay,az,sum;
  double w0,sumw,min_diff;
  double WMatomnum,longest_time;
  double ***List_T_in_Child;
  double **IMT,*weight;
  double *ko,*WMatomN;
  double xyz_c[4];

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /* set numprocs0 */
  if (numprocs<atomnum) numprocs0 = numprocs;
  else                  numprocs0 = atomnum;

  alloc_first_10 = 1;

  /* max_depth of level */

  if (numprocs==1 || atomnum==1)
    max_depth = 0;
  else  
    max_depth = (int)(log(((double)numprocs0-1.0+1.0e-7))/log(2.0)) + 1;

  /****************************
      allocation of arrays
  ****************************/

  Num_Procs_in_Child = (int**)malloc(sizeof(int*)*(max_depth+1));
  n = 1; 
  for (depth=0; depth<(max_depth+1); depth++){
    Num_Procs_in_Child[depth] = (int*)malloc(sizeof(int)*n);
    n *= 2;
  }

  Num_Atoms_in_Child = (int**)malloc(sizeof(int*)*(max_depth+1));
  n = 1; 
  for (depth=0; depth<(max_depth+1); depth++){
    Num_Atoms_in_Child[depth] = (int*)malloc(sizeof(int)*n);
    for (i=0; i<n; i++) Num_Atoms_in_Child[depth][i] = 0; 
    n *= 2;
  }

  IMT = (double**)malloc(sizeof(double*)*4);
  for (i=0; i<4; i++){
    IMT[i] = (double*)malloc(sizeof(double)*4);
  }

  ko = (double*)malloc(sizeof(double)*4);

  weight = (double*)malloc(sizeof(double)*(atomnum+1));

  /* set weight */    

  //if (weight_flag==0){
    for (i=1; i<=atomnum; i++){
      weight[i] = 1.0;
    }
  //}
  //else if (weight_flag==1){

  //  longest_time = 0.0;
  //  for (i=1; i<=atomnum; i++){
  //    if (longest_time<time_per_atom[i]) longest_time = time_per_atom[i];
  //  }

  //  for (i=1; i<=atomnum; i++){
  //    weight[i] = time_per_atom[i]/longest_time;
  //  }
  //}

  /* set Num_Procs_in_Child */    

  n = 2; 
  Num_Procs_in_Child[0][0] = numprocs0;

  for (depth=1; depth<(max_depth+1); depth++){
    for (i=0; i<n; i++){

      if (i%2==0){
        Num_Procs_in_Child[depth][i] = (Num_Procs_in_Child[depth-1][i/2]-1)/2+1;
      }  
      else{
        Num_Procs_in_Child[depth][i] = Num_Procs_in_Child[depth-1][i/2]-Num_Procs_in_Child[depth][i-1];  
      } 
    }
    n *= 2;
  }    

  /* set Num_Atoms_in_Child at depth=0 */    

  depth = 0; child = 0;
  Num_Atoms_in_Child[depth][child] = atomnum;

  /***************************************************************
   modified recursive bisection to set AN_in_Child at each depth 
  ***************************************************************/ 

  /**************************************************************************
   Since the size of the last index of List_AN_in_Child and List_T_in_Child
   is determined by the modified recursive bisection, they are allocated 
   on-the-fly. 
  **************************************************************************/

  /* allocation of List_AN_in_Child and List_T_in_Child */
  List_AN_in_Child = (int***)malloc(sizeof(int**)*(max_depth+1)); 
  List_T_in_Child = (double***)malloc(sizeof(double**)*(max_depth+1)); 
  List_AN_in_Child[0] = (int**)malloc(sizeof(int*)*1);
  List_T_in_Child[0] = (double**)malloc(sizeof(double*)*1);
  List_AN_in_Child[0][0] = (int*)malloc(sizeof(int)*(atomnum+1));
  List_T_in_Child[0][0]  = (double*)malloc(sizeof(double)*(atomnum+1));
  for (k=0; k<atomnum; k++)  List_AN_in_Child[depth][0][k] = k+1;

  n = 1; 
  for (depth=0; depth<max_depth; depth++){

    /* allocation of List_AN_in_Child and List_T_in_Child */
    List_AN_in_Child[depth+1] = (int**)malloc(sizeof(int*)*n*2);
    List_T_in_Child[depth+1] = (double**)malloc(sizeof(double*)*n*2);

    /**********************************************************
     reordering of atoms at depth using the inertia tensor
    **********************************************************/

    for (child=0; child<n; child++){

      /* get the number of atoms in the child */

      Na = Num_Atoms_in_Child[depth][child];

      /* calculate the centroid of atoms in the child */

      xyz_c[1] = 0.0; 
      xyz_c[2] = 0.0; 
      xyz_c[3] = 0.0; 

      for (k=0; k<Na; k++){
        m = List_AN_in_Child[depth][child][k];
        xyz_c[1] += Gxyz[m][1]*weight[m];
        xyz_c[2] += Gxyz[m][2]*weight[m];
        xyz_c[3] += Gxyz[m][3]*weight[m];
      }

      xyz_c[1] /= (double)Na;
      xyz_c[2] /= (double)Na;
      xyz_c[3] /= (double)Na;

      /* make inertia moment tensor */

      for (i=1; i<=3; i++){
        for (j=1; j<=3; j++){

          sum = 0.0;
          for (k=0; k<Na; k++){
	    m = List_AN_in_Child[depth][child][k];
	    sum += weight[m]*(Gxyz[m][i]-xyz_c[i])*(Gxyz[m][j]-xyz_c[j]);
	  }

          IMT[i][j] = sum;
        }
      }

      /* diagonalize the inertia moment tensor */

      Eigen_lapack(IMT,ko,3,3);

      /* find the principal axis */
  
      ax = IMT[1][3];
      ay = IMT[2][3];
      az = IMT[3][3];

      /* calculate the intervening variable, t */

      for (k=0; k<Na; k++){
        m = List_AN_in_Child[depth][child][k];
        t = ax*(Gxyz[m][1]-xyz_c[1]) + ay*(Gxyz[m][2]-xyz_c[2]) + az*(Gxyz[m][3]-xyz_c[3]);
        List_T_in_Child[depth][child][k] = t;
      }

      /* sorting atoms in the child based on t */

      qsort_double_int((long)Na,List_T_in_Child[depth][child],List_AN_in_Child[depth][child]);

      /* calculate the sum of weight in the child */

      sumw = 0.0;
      for (k=0; k<Na; k++){
        m = List_AN_in_Child[depth][child][k];
        sumw += weight[m];
      }

      /* find atomic index at which the bisection is made. */

      np = Num_Procs_in_Child[depth+1][2*child] + Num_Procs_in_Child[depth+1][2*child+1];
      w0 = (sumw*(double)Num_Procs_in_Child[depth+1][2*child])/(double)np;

      sumw = 0.0;
      min_diff = 10000000; 
      for (k=0; k<Na; k++){
        m = List_AN_in_Child[depth][child][k];
        sumw += weight[m];
        if (fabs(w0-sumw)<min_diff){
          min_diff = fabs(w0-sumw);
          k0 = k;
        }
      }

      /* adjust k0 to avoid the case that (# of atoms)<(# of processes) */

      if ( (((k0+1)<Num_Procs_in_Child[depth+1][2*child])
           ||
           ((Na-(k0+1))<Num_Procs_in_Child[depth+1][2*child+1]))
           && 
           1<Na ){

        k0 = Num_Procs_in_Child[depth+1][2*child] - 1; 
      }

      /* bisection of atoms in the child based on Num_Procs_in_Child */

      Num_Atoms_in_Child[depth+1][2*child  ] = k0 + 1;
      Num_Atoms_in_Child[depth+1][2*child+1] = Na - (k0+1);

      /* allocation of List_AN_in_Child and List_T_in_Child */
      List_AN_in_Child[depth+1][2*child]   = (int*)malloc(sizeof(int)*Num_Atoms_in_Child[depth+1][2*child]);
      List_T_in_Child[depth+1][2*child]   = (double*)malloc(sizeof(double)*Num_Atoms_in_Child[depth+1][2*child]);
      List_AN_in_Child[depth+1][2*child+1] = (int*)malloc(sizeof(int)*Num_Atoms_in_Child[depth+1][2*child+1]);
      List_T_in_Child[depth+1][2*child+1] = (double*)malloc(sizeof(double)*Num_Atoms_in_Child[depth+1][2*child+1]);

      /* copy depth -> depth+1 */

      for (k=0; k<Num_Atoms_in_Child[depth+1][2*child]; k++){
        List_AN_in_Child[depth+1][2*child][k] = List_AN_in_Child[depth][child][k];
      }

      for (k=0; k<Num_Atoms_in_Child[depth+1][2*child+1]; k++){
        m = Num_Atoms_in_Child[depth+1][2*child]+k;
        List_AN_in_Child[depth+1][2*child+1][k] = List_AN_in_Child[depth][child][m];
      }

    } /* child */     

    /* doubling of n */
    n *= 2;

  } /* depth */

  /*
  if (myid==0){
    n = 1; 
    for (depth=0; depth<=max_depth; depth++){
      for (child=0; child<n; child++){

	Na = Num_Atoms_in_Child[depth][child];

	for (k=0; k<Na; k++){
	  m = List_AN_in_Child[depth][child][k];
	  t = List_T_in_Child[depth][child][k];
	  printf("depth=%2d child=%2d k=%2d m=%2d t=%15.12f\n",depth,child,k,m,t);
	}
      }
      n *= 2;
    }
  }
  */

  /***************************************************************
                 allocation of atoms to processes
  ***************************************************************/ 

  /*
  sorting of atoms in each child at max_depth. 
  if the sorting is not performed, the force calculations
  related to HNL3 and 4B will be failed.
  */

  n = 1; 
  for (depth=0; depth<max_depth; depth++) n *= 2;

  for (child=0; child<n; child++){
    Na = Num_Atoms_in_Child[max_depth][child];
    qsort_int1((long)Na,List_AN_in_Child[depth][child]);
  }  

  /* set G2ID, M2G, Matomnum, and WMatomnum */

  n = 1; 
  for (depth=0; depth<max_depth; depth++) n *= 2;

  Matomnum = 0;
  WMatomnum = 0.0;
  ID = 0;

  for (child=0; child<n; child++){

    Na = Num_Atoms_in_Child[max_depth][child];

    if (Na!=0){ 

      for (k=0; k<Na; k++){
        m = List_AN_in_Child[max_depth][child][k];
        G2ID[m] = ID; 
      }

      if (myid==ID){

        Matomnum = Na; 
        if (alloc_first_10==0) free(M2G);
	M2G = (int*)malloc(sizeof(int)*(Matomnum+2));
	alloc_first_10 = 0;

        for (k=0; k<Na; k++){
          m = List_AN_in_Child[max_depth][child][k];
          M2G[k+1] = m;
        }

        WMatomnum = 0.0;
        for (k=0; k<Na; k++){
          m = List_AN_in_Child[max_depth][child][k];
          WMatomnum += weight[m];
        }
      }

      ID++; 
    }
  }

  /****************************************
    find Max_Matomnum, MatomN and WMatomN
  ****************************************/

  MatomN = (int*)malloc(sizeof(int)*numprocs);
  WMatomN = (double*)malloc(sizeof(double)*numprocs);

  MatomN[myid]  = Matomnum;
  WMatomN[myid] = WMatomnum;

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&MatomN[ID],  1, MPI_INT, ID, mpi_comm_level1);
    MPI_Bcast(&WMatomN[ID], 1, MPI_DOUBLE, ID, mpi_comm_level1);
  } 

  /* find Max_Matomnum */

  //Max_Matomnum = 0;
  //for (ID=0; ID<numprocs; ID++){
  //  if (Max_Matomnum<MatomN[ID]) Max_Matomnum = MatomN[ID];
  //}     

  /*********************************************
              print the information 
  *********************************************/

  if (myid==Host_ID && 0<level_stdout){
    printf("\n");
    printf("*******************************************************\n"); 
    printf("  Allocation of atoms to proccesors \n");
    printf("*******************************************************\n\n"); 
  }

  for (ID=0; ID<numprocs; ID++){

    if (myid==Host_ID && 0<level_stdout){
      printf(" proc = %3d  # of atoms=%4d  estimated weight=%16.5f\n",
	     ID,MatomN[ID],WMatomN[ID]);
    }
  }     

  if (myid==Host_ID && 0<level_stdout) printf("\n\n\n");

  /****************************************
         initialize time_per_atom
  ****************************************/

  //for (k=1; k<=atomnum; k++) time_per_atom[k] = 0.0;

  /*
  if (myid==Host_ID){

    int i,k;
    char AtomName[38][10]=
      {"E","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
       "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
       "Ga","Ge","As","Se","Br","Kr"
      };

    for (i=1; i<=atomnum; i++){
      k = (G2ID[i])%36+4;
      printf("%5s  %15.12f %15.12f %15.12f  0.0 0.0 0.0\n",
               AtomName[k],Gxyz[i][1]*BohrR,Gxyz[i][2]*BohrR,Gxyz[i][3]*BohrR);
    }
  }
  MPI_Finalize();
  exit(0);
  */
  
  /*
  {
    int i1,i2,k1,k2;
    double dx,dy,dz,r,rmax;

    rmax = 0.0;
    for (i1=1; i1<=Matomnum; i1++){
      k1 = M2G[i1];

      for (i2=1; i2<=Matomnum; i2++){
        k2 = M2G[i2];

        dx = Gxyz[k1][1] - Gxyz[k2][1];
        dy = Gxyz[k1][2] - Gxyz[k2][2];
        dz = Gxyz[k1][3] - Gxyz[k2][3];
        r = sqrt(dx*dx+dy*dy+dz*dz);
        if (rmax<r) rmax = r;
      }
    }

    printf("ABC1 myid=%2d rmax=%15.12f\n",myid,rmax);fflush(stdout);
  }

  MPI_Finalize();
  exit(0);
  */


  /****************************************
            freeing of arrays
  ****************************************/

  free(WMatomN);
  free(MatomN);

  free(List_AN_in_Child[0][0]);
  free(List_T_in_Child[0][0]);
  free(List_AN_in_Child[0]);
  free(List_T_in_Child[0]);

  n = 1; 
  for (depth=0; depth<max_depth; depth++){

    for (i=0; i<n*2; i++){
      free(List_T_in_Child[depth+1][i]);
      free(List_AN_in_Child[depth+1][i]);
    }
    free(List_T_in_Child[depth+1]);
    free(List_AN_in_Child[depth+1]);

    n *= 2;
  }
  free(List_T_in_Child);
  free(List_AN_in_Child);

  free(weight);
  free(ko);

  for (i=0; i<4; i++){
    free(IMT[i]);
  }
  free(IMT);

  for (depth=0; depth<(max_depth+1); depth++){
    free(Num_Atoms_in_Child[depth]);
  }
  free(Num_Atoms_in_Child);

  for (depth=0; depth<(max_depth+1); depth++){
    free(Num_Procs_in_Child[depth]);
  }
  free(Num_Procs_in_Child);
}

   
void Eigen_lapack(double **a, double *ko, int n, int EVmax)
{
  static int solver_flag=0;
  int i,j,k,po,okay_flag,iterN;
  double sum;
  double **b;

  b=(double**)malloc(sizeof(double*)*(n+1));
  for (i=0; i<(n+1); i++){
    b[i]=(double*)malloc(sizeof(double)*(n+1));
  }

  for (i=0; i<(n+1); i++){
    for (j=0; j<(n+1); j++){
      b[i][j] = a[i][j]; 
    }
  }

  iterN = 1;
  okay_flag = 0;

  do {

    for (i=0; i<(n+1); i++){
      for (j=0; j<(n+1); j++){
	a[i][j] = b[i][j]; 
      }
    }

    if      (solver_flag==0) Eigen_lapack_x(a, ko, n, EVmax);
    else if (solver_flag==1) Eigen_HH(a, ko, n, EVmax);  
    else if (solver_flag==2) Eigen_lapack_d(a, ko, n, EVmax); 

    po = 0; 

    i = 1;
    do {
      j = 1;
      do {

	sum = 0.0;
	for (k=1; k<=n; k++){
	  sum += a[i][k]*a[j][k]; 
	}

	if      (i==j && 0.00001<fabs(sum-1.0)){ po = 1; }
	else if (i!=j && 0.00001<fabs(sum))    { po = 1; }

        j = j + 20;

      } while (po==0 && j<=EVmax); 

      i = i + 20;

    } while (po==0 && i<=EVmax); 

    if (po==1){
      solver_flag++; 
      solver_flag = solver_flag % 3;
    }
    else {
      okay_flag = 1;
    }

    /*
    printf("iterN=%2d solver_flag=%2d po=%2d okay_flag=%2d\n",iterN,solver_flag,po,okay_flag);
    */

    iterN++;

  } while (okay_flag==0 && iterN<4);  

  for (i=0; i<(n+1); i++){
    free(b[i]);
  }
  free(b);

} 



void Eigen_lapack_x(double **a, double *ko, int n0, int EVmax)
{

  /*
    F77_NAME(dsyevx,DSYEVX)()
  
    input:  n;
    input:  a[n][n];  matrix A
    output: a[n][n];  eigevectors
    output: ko[n];    eigenvalues 
  */
    
  char *name="Eigen_lapack_x";

  char  *JOBZ="V";
  char  *RANGE="I";
  char  *UPLO="L";

  INTEGER n=n0;
  INTEGER LDA=n0;
  double VL,VU; /* dummy */
  INTEGER IL,IU; 
  double ABSTOL=LAPACK_ABSTOL;
  INTEGER M;

  double *A,*Z;
  INTEGER LDZ=n;
  INTEGER LWORK;
  double *WORK;
  INTEGER *IWORK;
  INTEGER *IFAIL, INFO;

  int i,j;

  A=(double*)malloc(sizeof(double)*n*n);
  Z=(double*)malloc(sizeof(double)*n*n);

  LWORK=n*8;
  WORK=(double*)malloc(sizeof(double)*LWORK);
  IWORK=(INTEGER*)malloc(sizeof(INTEGER)*n*5);
  IFAIL=(INTEGER*)malloc(sizeof(INTEGER)*n);

  IL = 1;
  IU = EVmax;
 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
       A[i*n+j] = a[i+1][j+1];
    }
  }

#if 0
  printf("A=\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       printf("%f ",A[i*n+j]);
    }
    printf("\n");
  }
  fflush(stdout);
#endif

  F77_NAME(dsyevx,DSYEVX)( JOBZ, RANGE, UPLO, &n, A, &LDA, &VL, &VU, &IL, &IU,
           &ABSTOL, &M, ko, Z, &LDZ, WORK, &LWORK, IWORK,
           IFAIL, &INFO ); 


  /* store eigenvectors */
  for (i=0;i<EVmax;i++) {
    for (j=0;j<n;j++) {
      /*  a[i+1][j+1]= Z[i*n+j]; */
      a[j+1][i+1]= Z[i*n+j];
    }
  }

  /* shift ko by 1 */
  for (i=EVmax; i>=1; i--){
    ko[i]= ko[i-1];
  }

  /*
  if (INFO>0) {
    printf("\n%s: error in dsyevx_, info=%d\n\n",name,INFO);
  }
  */

  if (INFO<0) {
     printf("%s: info=%d\n",name,INFO);
     exit(10);
  }
   
  free(IFAIL); free(IWORK); free(WORK); free(Z); free(A);

}


void Eigen_lapack_d(double **a, double *ko, int n0, int EVmax)
{

  /* 
    F77_NAME(dsyevd,DSYEVD)()
  
    input:  n;
    input:  a[n][n];  matrix A
    output: a[n][n];  eigevectors
    output: ko[n];    eigenvalues 
  */
    
  static char *name="Eigen_lapack_d";

  char  *JOBZ="V";
  char  *UPLO="L";

  INTEGER n=n0;
  INTEGER LDA=n;
  double VL,VU; /* dummy */
  INTEGER IL,IU; 
  double ABSTOL=LAPACK_ABSTOL;
  INTEGER M;

  double *A;
  INTEGER LDZ=n;
  INTEGER LWORK,LIWORK;
  double *WORK;
  INTEGER *IWORK;
  INTEGER INFO;

  int i,j;

  A=(double*)malloc(sizeof(double)*n*n);

  LWORK=  1 + 6*n + 2*n*n;
  WORK=(double*)malloc(sizeof(double)*LWORK);

  LIWORK = 3 + 5*n;
  IWORK=(INTEGER*)malloc(sizeof(INTEGER)*LIWORK);


  IL = 1;
  IU = EVmax; 
 
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       A[i*n+j]= a[i+1][j+1];
    }
  }

#if 0
  printf("A=\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       printf("%f ",A[i*n+j]);
    }
    printf("\n");
  }
  fflush(stdout);
#endif

  F77_NAME(dsyevd,DSYEVD)( JOBZ, UPLO, &n, A, &LDA, ko, WORK, &LWORK, IWORK, &LIWORK, &INFO ); 

  /* store eigenvectors */
  for (i=0;i<EVmax;i++) {
    for (j=0;j<n;j++) {
      /* a[i+1][j+1]= Z[i*n+j]; */
      a[j+1][i+1]= A[i*n+j];
    }
  }

  /* shift ko by 1 */
  for (i=EVmax; i>=1; i--){
    ko[i]= ko[i-1];
  }

  /*
  if (INFO>0) {
     printf("\n%s: error in dsyevd_, info=%d\n\n",name,INFO);
  }
  */

  if (INFO<0) {
     printf("%s: info=%d\n",name,INFO);
     exit(10);
  }
   
  free(IWORK); free(WORK); free(A);
}



#pragma optimization_level 1
void Eigen_HH(double **ac, double *ko, int n, int EVmax)
{
  /**********************************************************************
    Eigen_HH:

    Eigen_HH.c is a subroutine to solve a seqular equation without an
    overlap matrix using Householder method and lapack's dstevx, dstegr, 
    or dstedc.

    Log of Eigen_HH.c:

       Nov/22/2004  Released by T.Ozaki

  ***********************************************************************/

  double ABSTOL=LAPACK_ABSTOL;
  double **ad,*D,*E,
                *b1,*u,*uu,
                *p,*q,*s,*c,
                s1,s2,s3,ss,u1,u2,r,p1,
                s20,s21,s22,s23,
                xsum,bunbo,si,co,sum,
                a1,a2,a3,a4,a5,a6,b7,
                x1,x2,xap,tmp1,tmp2,
                bb,bb1,ui,uj,uij;

  double ss0,ss1,ss2,ss3;
  double r0,r1,r2,r3,p10,p11,p12,p13;
  double tmp10,tmp11,tmp12,tmp13;
  double tmp20,tmp21,tmp22,tmp23;
 
  int jj,jj1,jj2,k,ii,ll,i3,i2,j2,i1s,
             i,j,i1,j1,n1,n2,ik,
             jk,po1,nn,count,ks;

  double Stime, Etime;
  double Stime1, Etime1;
  double Stime2, Etime2;
  double time1,time2;

  //int measure_time;

  //measure_time=0;

  /****************************************************
    allocation of arrays:
  ****************************************************/

  n2 = n + 5;

  ad = (double**)malloc(sizeof(double*)*n2);
  for (i=0; i<n2; i++){
    ad[i] = (double*)malloc(sizeof(double)*n2);
  }

  b1 = (double*)malloc(sizeof(double)*n2);
  u = (double*)malloc(sizeof(double)*n2);
  uu = (double*)malloc(sizeof(double)*n2);
  p = (double*)malloc(sizeof(double)*n2);
  q = (double*)malloc(sizeof(double)*n2);
  s = (double*)malloc(sizeof(double)*n2);
  c = (double*)malloc(sizeof(double)*n2);

  D = (double*)malloc(sizeof(double)*n2);
  E = (double*)malloc(sizeof(double)*n2);

  for (i=1; i<=(n+2); i++){
    uu[i] = 0.0;
  }

  //if (measure_time==1) printf("size n=%3d EVmax=%2d\n",n,EVmax);
  //if (measure_time==1) dtime(&Stime);

  /****************************************************
                   Householder method
  ****************************************************/

  for (i=1; i<=(n-2); i++){

    /* original version */

    /*
    s1 = ac[i+1][i] * ac[i+1][i];
    s2 = 0.0;
    u[i+1] = ac[i+1][i];
    for (i1=i+2; i1<=n; i1++){
      tmp1 = ac[i1][i]; 
      s2 += tmp1*tmp1;
      u[i1] = tmp1;
    }
    s3 = fabs(s1 + s2);
    */

    /* unrolling version */

    s1 = ac[i+1][i] * ac[i+1][i];
    u[i+1] = ac[i+1][i];

    s20 = 0.0;
    s21 = 0.0;
    s22 = 0.0;
    s23 = 0.0;

    for (i1=i+2; i1<=(n-3); i1+=4){

      s20 += ac[i][i1+0]*ac[i][i1+0];
      s21 += ac[i][i1+1]*ac[i][i1+1];
      s22 += ac[i][i1+2]*ac[i][i1+2];
      s23 += ac[i][i1+3]*ac[i][i1+3];

      u[i1+0] = ac[i][i1+0];
      u[i1+1] = ac[i][i1+1];
      u[i1+2] = ac[i][i1+2];
      u[i1+3] = ac[i][i1+3];
    }

    i1s = n + 1 - (n+1-(i+2))%4;

    for (i1=i1s; ((i+2)<=i1 && i1<=n); i1++){
      tmp1 = ac[i][i1]; 
      s20 += tmp1*tmp1;
      u[i1] = tmp1;
    }

    s2 = s20 + s21 + s22 + s23;
    s3 = fabs(s1 + s2);

    if (ABSTOL<fabs(ac[i+1][i])){
      if (ac[i+1][i]<0.0)    s3 =  sqrt(s3);
      else                   s3 = -sqrt(s3);
    }
    else{
      s3 = sqrt(s3);
    }

    if (ABSTOL<fabs(s2)){

      ss = ac[i+1][i];
      ac[i+1][i] = s3;
      ac[i][i+1] = s3;
      u[i+1] = u[i+1] - s3;
      u1 = s3 * s3 - ss * s3;
      u2 = 2.0 * u1;
      uu[i] = u2;
      b1[i] = ss - s3;

      /* original version */

      /*
      r = 0.0;
      for (i1=i+1; i1<=n; i1++){
	p1 = 0.0;
	for (j=i+1; j<=n; j++){
	  p1 += ac[i1][j] * u[j];
	}
	p[i1] = p1 / u1;
	r += u[i1] * p[i1];
      }
      r = r / u2;
      */

      /* unrolling version */

      r0 = 0.0;
      r1 = 0.0;
      r2 = 0.0;
      r3 = 0.0;

      for (i1=i+1; i1<=(n-3); i1+=4){

	p10 = 0.0;
	p11 = 0.0;
	p12 = 0.0;
	p13 = 0.0;

	for (j=i+1; j<=n; j++){
	  p10 += ac[i1+0][j] * u[j];
	  p11 += ac[i1+1][j] * u[j];
	  p12 += ac[i1+2][j] * u[j];
	  p13 += ac[i1+3][j] * u[j];
	}

	p[i1+0] = p10 / u1;
	p[i1+1] = p11 / u1;
	p[i1+2] = p12 / u1;
	p[i1+3] = p13 / u1;

	r0 += u[i1+0] * p[i1+0];
	r1 += u[i1+1] * p[i1+1];
	r2 += u[i1+2] * p[i1+2];
	r3 += u[i1+3] * p[i1+3];
      }

      i1s = n + 1 - (n+1-(i+1))%4;

      for (i1=i1s; ((i+1)<=i1 && i1<=n); i1++){
	p1 = 0.0;
	for (j=i+1; j<=n; j++){
	  p1 += ac[i1][j] * u[j];
	}
	p[i1] = p1 / u1;
	r0 += u[i1] * p[i1];
      }

      r = (r0+r1+r2+r3) / u2;

      /* original version */

      /*
      for (i1=i+1; i1<=n; i1++){
	q[i1] = p[i1] - r * u[i1];
      }
      */

      /* unrolling version */

      for (i1=i+1; i1<=(n-3); i1+=4){
	q[i1+0] = p[i1+0] - r * u[i1+0];
	q[i1+1] = p[i1+1] - r * u[i1+1];
	q[i1+2] = p[i1+2] - r * u[i1+2];
	q[i1+3] = p[i1+3] - r * u[i1+3];
      }

      i1s = n + 1 - (n+1-(i+1))%4;

      for (i1=i1s; ((i+1)<=i1 && i1<=n); i1++){
	q[i1] = p[i1] - r * u[i1];
      }

      /* original version */

      /*
      for (i1=i+1; i1<=n; i1++){
        tmp1 = u[i1];
        tmp2 = q[i1]; 
	for (j1=i+1; j1<=n; j1++){
	  ac[i1][j1] -= tmp1 * q[j1] + tmp2 * u[j1];
	}
      }
      */

      /* unrolling version */

      for (i1=i+1; i1<=(n-3); i1+=4){

        tmp10 = u[i1+0];
        tmp11 = u[i1+1];
        tmp12 = u[i1+2];
        tmp13 = u[i1+3];

        tmp20 = q[i1+0]; 
        tmp21 = q[i1+1]; 
        tmp22 = q[i1+2]; 
        tmp23 = q[i1+3]; 

	for (j1=i+1; j1<=n; j1++){
	  ac[i1+0][j1] -= tmp10 * q[j1] + tmp20 * u[j1];
	  ac[i1+1][j1] -= tmp11 * q[j1] + tmp21 * u[j1];
	  ac[i1+2][j1] -= tmp12 * q[j1] + tmp22 * u[j1];
	  ac[i1+3][j1] -= tmp13 * q[j1] + tmp23 * u[j1];
	}
      }

      i1s = n + 1 - (n+1-(i+1))%4;

      for (i1=i1s; ((i+1)<=i1 && i1<=n); i1++){
        tmp1 = u[i1];
        tmp2 = q[i1]; 
	for (j1=i+1; j1<=n; j1++){
	  ac[i1][j1] -= tmp1 * q[j1] + tmp2 * u[j1];
	}
      }

    }
  }

  for (i=1; i<=n; i++){
    for (j=1; j<=n; j++){
      ad[i][j] = ac[i][j];
    }
  }

  /*
  for (i=1; i<=n; i++){
    printf("i=%4d  D=%18.15f E=%18.15f\n",i,ad[i][i],ad[i][i+1]);
  }
  */

  //if (measure_time==1){
  //  dtime(&Etime);
  //  printf("T1   %15.12f\n",Etime-Stime);
  //}

  /****************************************************
                  call a lapack routine
  ****************************************************/

  //if (measure_time==1) dtime(&Stime);

  for (i=1; i<=n; i++){
    D[i-1] = ad[i][i];
    E[i-1] = ad[i][i+1];
  }

  lapack_dstegr1(n,EVmax,D,E,ko,ac);

  //if      (dste_flag==0) lapack_dstegr1(n,EVmax,D,E,ko,ac);
  //else if (dste_flag==1) lapack_dstedc1(n,D,E,ko,ac);
  //else if (dste_flag==2) lapack_dstevx1(n,EVmax,D,E,ko,ac);

  //if (measure_time==1){
  //  dtime(&Etime);
  //  printf("T2   %15.12f\n",Etime-Stime);
  //}

  /****************************************************
    transformation of eigenvectors to original space
  ****************************************************/

  //if (measure_time==1) dtime(&Stime);

  for (i=2; i<=n-1; i++){
    ad[i-1][i] = b1[i-1];
  }

  /* original version */

  /*
  for (k=1; k<=EVmax; k++){
    for (nn=2; nn<=n-1; nn++){
      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){
	ss = 0.0;
	for (i=n-nn+1; i<=n; i++){
	  ss += ad[n-nn][i] * ac[k][i];
	}
	ss = 2.0*ss/uu[n-nn];
	for (i=n-nn+1; i<=n; i++){
	  ac[k][i] -= ss * ad[n-nn][i];
	}
      }
    }
  }
  */

  /* unrolling version */

  for (k=1; k<=(EVmax-3); k+=4){

    for (nn=2; nn<=n-1; nn++){
      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){

	ss0 = 0.0;
	ss1 = 0.0;
	ss2 = 0.0;
	ss3 = 0.0;

	for (i=n-nn+1; i<=n; i++){
	  ss0 += ad[n-nn][i] * ac[k+0][i];
	  ss1 += ad[n-nn][i] * ac[k+1][i];
	  ss2 += ad[n-nn][i] * ac[k+2][i];
	  ss3 += ad[n-nn][i] * ac[k+3][i];
	}

	ss0 = 2.0*ss0/uu[n-nn];
	ss1 = 2.0*ss1/uu[n-nn];
	ss2 = 2.0*ss2/uu[n-nn];
	ss3 = 2.0*ss3/uu[n-nn];

	for (i=n-nn+1; i<=n; i++){
	  ac[k+0][i] -= ss0 * ad[n-nn][i];
	  ac[k+1][i] -= ss1 * ad[n-nn][i];
	  ac[k+2][i] -= ss2 * ad[n-nn][i];
	  ac[k+3][i] -= ss3 * ad[n-nn][i];
	}

      }
    }
  }

  ks = EVmax - EVmax%4 + 1;

  for (k=ks; k<=EVmax; k++){
    for (nn=2; nn<=n-1; nn++){
      if ( (1.0e-3*ABSTOL)<fabs(uu[n-nn])){
	ss = 0.0;
	for (i=n-nn+1; i<=n; i++){
	  ss += ad[n-nn][i] * ac[k][i];
	}
	ss = 2.0*ss/uu[n-nn];
	for (i=n-nn+1; i<=n; i++){
	  ac[k][i] -= ss * ad[n-nn][i];
	}
      }
    }
  }

  //if (measure_time==1){
  //  dtime(&Etime);
  //  printf("T4   %15.12f\n",Etime-Stime);
  //}

  /****************************************************
                     normalization
  ****************************************************/

  //if (measure_time==1) dtime(&Stime);

  for (j=1; j<=EVmax; j++){
    sum = 0.0;
    for (i=1; i<=n; i++){
      sum = sum + ac[j][i] * ac[j][i];
    }
    sum = 1.0/sqrt(sum);
    for (i=1; i<=n; i++){
      ac[j][i] = ac[j][i] * sum;
    }
  }

  //if (measure_time==1){
  //  dtime(&Etime);
  //  printf("T5   %15.12f\n",Etime-Stime);
  //}

  /****************************************************
            Eigenvectors to the "ac" array
  ****************************************************/

  for (i=1; i<=n; i++){
    for (j=(i+1); j<=n; j++){
      tmp1 = ac[i][j];
      tmp2 = ac[j][i];
      ac[i][j] = tmp2;
      ac[j][i] = tmp1;
    }
  }

  /****************************************************
                  freeing of arrays:
  ****************************************************/

  for (i=0; i<n2; i++){
    free(ad[i]);
  }
  free(ad);

  free(b1);
  free(u);
  free(uu);
  free(p);
  free(q);
  free(s);
  free(c);
  free(D);
  free(E);
}

