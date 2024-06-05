   
#define PI              3.1415926535897932384626
#define BYTESIZE        8                        /* Don't change!! */
#define kB              0.00008617251324000000   /* eV/K           */          
#define BohrR           0.529177249              /* Angstrom       */
#define eV2Hartree      27.2113845                

#define Host_ID             0        /* ID of the host CPU in MPI */
#define YOUSO10     500        /* length of a file name                           */
#define NYOUSO       60        /* # of YOUSO                                      */
#define YOUSO26      61         /* size of the second index in Gxyz[][]           */
#define Supported_MaxL      4        /* supported max angular momentum for basis orbital */
#define LAPACK_ABSTOL     6.0e-15    /* absolute error tolerance for lapack routines */
typedef double     Type_Orbs_Grid;       /* type of Orbs_Grid */

#ifndef __sqr_definition___
#define sqr(x)   ( (x)*(x) )
#define __sqr_definition___
#endif

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif

#ifdef nompi

#ifndef ___MPI_Comm_definition___
typedef int MPI_Comm;
#define ___MPI_Comm_definition___ 
#endif

#ifndef ___MPI_Status_definition___
typedef struct MPIStatus{int i;}  MPI_Status;  
#define ___MPI_Status_definition___ 
#endif

#ifndef ___MPI_Request_definition___
typedef struct MPIRequest{int i;} MPI_Request;  
#define ___MPI_Request_definition___ 
#endif

#else
#include "mpi.h"
#endif

MPI_Comm  mpi_comm_level1;
MPI_Comm  MPI_COMM_WORLD1;
int Num_Procs;
int NUMPROCS_MPI_COMM_WORLD,MYID_MPI_COMM_WORLD;

/****************************/
/* input data from flpq.inp */
/****************************/
double scf_fixed_dens_origin[3]; /* fukuda */
int Dens_Ngrid[3];
char DirnameLPQ[100];
int flag_DMkmu,flag_DMmu;
int flag_DMk;
int flag_DMmu_range;
int flag_grid_vec;
int flag_save_memory_3d;
int loop_Ngrid1, N_loop_Ngrid1;
int Num_DM_orb,Num_DM_k,Num_DM_k_orb,Num_DM_spin;


/**********************************/
/* input data from OpenMX_LPQ.dat */
/**********************************/
int SpinP_switch;
int Solver;
int Cnt_switch;
int Cnt_kind;
int Calc_CntOrbital_ON;
int MD_switch;
int level_stdout;
int ESM_switch;
int TCpyCell,CpyCell;

int Ngrid1, Ngrid2, Ngrid3;

int List_YOUSO[NYOUSO];
//List_YOUSO[7]
//List_YOUSO[11] /* = (int)(Max_GridN_Atom*ScaleSize) + 1; */
//List_YOUSO[12] /*  = (int)(Max_NumOLG*ScaleSize) + 1;    */
//List_YOUSO[17] /* = (int)(Max_OneD_Grids*ScaleSize);     */
//List_YOUSO[18]
//List_YOUSO[21]
//List_YOUSO[24]
//List_YOUSO[25]


int SpeciesNum;
int Max_FSNAN;
int Max_NumOLG;
double ScaleSize;

int Matomnum, MatomnumF, atomnum;

int NN_A2B_S,NN_A2B_R;

double tv[4][4],rtv[4][4];
double grid_vec[4][4];
double tv_ori[4][4];

int flag_energy_range_DM;
double DM_energy_range[2];



/*******************************************************
 int *Spe_Total_NO; 
 the number of primitive atomic orbitals in a species
  size: Spe_Total_NO[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *Spe_Total_NO;

/*******************************************************
 int *Spe_Total_CNO; 
 the number of contracted atomic orbitals in a species
  size: Spe_Total_CNO[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *Spe_Total_CNO;


/*******************************************************
 int *Spe_Num_Mesh_PAO; 
 the number of grids for atomic orbitals
  size: Spe_Num_Mesh_PAO[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *Spe_Num_Mesh_PAO;

/*******************************************************
 double **Spe_PAO_RV;
  logarithmic radial mesh (r=exp(r)) for PAO 
  size: Spe_PAO_XV[List_YOUSO[18]]
                  [List_YOUSO[21]]
  allocation: call as Allocate_Arrays(6) in SetPara_DFT.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **Spe_PAO_RV;

/*******************************************************
 double ****Spe_PAO_RWF;
  radial parts of basis orbitals on radial mesh of PAO 
  size: Spe_PAO_RWF[List_YOUSO[18]]
                   [List_YOUSO[25]+1]
                   [List_YOUSO[24]]
                   [List_YOUSO[21]]
  allocation: call as Allocate_Arrays(6) in SetPara_DFT.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ****Spe_PAO_RWF;

/*******************************************************
 int *Spe_MaxL_Basis; 
 the maximum "l" component of used atomic orbitals inv
 each species
  size: Spe_MaxL_Basis[SpeciesNum]
*******************************************************/
int *Spe_MaxL_Basis;

/*******************************************************
 int **Spe_Num_Basis; 
 the number of multiplicity of primitive radial parts
 for each "l" component in an species
  size: Spe_Num_Basis[SpeciesNum][6]
  allocation: call as Allocate_Arrays(0) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Spe_Num_Basis;

/*******************************************************
 int **Spe_Num_CBasis; 
 the number of multiplicity of contracted radial parts
 for each "l" component in an species
  size: Spe_Num_CBasis[SpeciesNum][6]
  allocation: call as Allocate_Arrays(0) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Spe_Num_CBasis;

/*******************************************************
 int **Spe_Specified_Num;
  a table which converts index of contracted orbitals
  to that of primitive orbitals
  size: Spe_Specified_Num[List_YOUSO[18]]
                         [Spe_Total_NO[spe]]  
  allocation: in Set_BasisPara() of SetPara_DFT.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Spe_Specified_Num;

/*******************************************************
 int ***Spe_Trans_Orbital;
  a table which converts index of contracted orbitals
  to that of primitive orbitals
  size: Spe_Trans_Orbital[List_YOUSO[18]]
                         [Spe_Total_NO[spe]]  
                         [List_YOUSO[24]]
  allocation: in Set_BasisPara() of SetPara_DFT.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int ***Spe_Trans_Orbital;

/*******************************************************
 double *Spe_Atom_Cut1; 
 cutoff radius of atomic orbitals for each species
  size: Spe_Atom_Cut1[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double *Spe_Atom_Cut1;

/*******************************************************
 int *Spe_WhatAtom; 
 atomic number in the periodic table for each species
  size: Spe_WhatAtom[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *Spe_WhatAtom;

/*******************************************************
 double *InitN_USpin; 
  the number of the upspin electon of initial atoms
  size: InitN_USpin[atomnum+1]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double *InitN_USpin;

/*******************************************************
 double *InitN_DSpin; 
  the number of the upspin electon of initial atoms
  size: InitN_DSpin[atomnum+1]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double *InitN_DSpin;

/*******************************************************
 double *Spe_Core_Charge; 
 effective core charge of each species
  size: Spe_Core_Charge[SpeciesNum]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double *Spe_Core_Charge;

/*******************************************************
 int *M2G
  M2G gives a conversion from the medium
  index to the global indicies of atoms.
  size: M2G[Matomnum+1]
  allocation: in Set_Allocate_Atom2CPU.c
  free:       call as Free_Arrays(0) in openmx.c
              and in Set_Allocate_Atom2CPU.c
*******************************************************/
int *M2G;

/*******************************************************
 int *F_M2G; 
  F_M2G gives a conversion from the medium
  index (Matomnum+MatomnumF)
  to the global indicies of atoms.
  size: F_M2G[Matomnum+MatomnumF+1],
  allocation: in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
              and in truncation.c
*******************************************************/
int *F_M2G;

/*******************************************************
 int *F_G2M;

  F_G2M give a conversion from the
  global atom number to the medium atom number
  for atoms sent from ID in the size of
  F_Rcv_Num[ID]. 
  size: F_G2M[atomnum+1]
  allocation: Allocation_Arrays(1) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *F_G2M;

/*******************************************************
 int *WhatSpecies; 
 array to specify species for each atom in the system 
  size: WhatSpecies[atomnum+1]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *WhatSpecies;

/*******************************************************
 int *GridN_Atom; 
 the number of grids overlaping to each atom
  size: GridN_Atom[atomnum+1]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *GridN_Atom;

/*******************************************************
 int *FNAN; 
 the number of first neighboring atoms
  size: FNAN[atomnum+1]
  allocation: call as Allocate_Arrays(2) in readfile.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *FNAN;

/*******************************************************
 int **ncn; 
  grobal index number for cell of neighboring atoms of
  an atom ct_AN
  size: ncn[atomnum+1][Max_FSNAN*ScaleSize+1]
  allocation: call as Allocate_Arrays(3) in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **ncn;

/*******************************************************
 int *G2ID;
  G2ID gives a proccesor ID allocated to each atom
  with a global atom index.
  size: G2ID[atomnum+1];
  allocation: Allocation_Arrays(1) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *G2ID;

/*******************************************************
 int *F_Snd_Num;

  F_Snd_Num gives the number of atoms of which informations,
  related by FNAN, are transfered from myid to ID.
  size: F_Snd_Num[numprocs]
  allocation: Allocation_Arrays(0) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *F_Snd_Num;

/*******************************************************
 int *S_Snd_Num;

  S_Snd_Num gives the number of atoms of which informations,
  related by SNAN, are transfered from myid to ID.
  size: S_Snd_Num[numprocs]
  allocation: Allocation_Arrays(0) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *S_Snd_Num;

/*******************************************************
 int *F_Rcv_Num;

  F_Rcv_Num gives the number of atoms of which informations,
  related by FNAN, are recieved at myid from ID.
  size: F_Rcv_Num[numprocs]
  allocation: Allocation_Arrays(0) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *F_Rcv_Num;

/*******************************************************
 int *S_Rcv_Num;

  S_Rcv_Num gives the number of atoms of which informations,
  related by SNAN, are recieved at myid from ID.
  size: S_Rcv_Num[numprocs]
  allocation: Allocation_Arrays(0) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *S_Rcv_Num;

/*******************************************************
 int *F_TopMAN;

  F_TopMAN and S_TopMAN give the first medium
  atom number in atoms sent from ID in the size of
  F_Rcv_Num[ID],
  respectively.
  size: F_TopMAN[numprocs]
  allocation: Allocation_Arrays(0) in Input_std()
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int *F_TopMAN;

/*******************************************************
 int **Snd_MAN;
  Snd_MAN is a medium atom index of which informations
  are sent to a processor ID.
  size: Snd_MAN[numprocs][FS_Snd_Num[ID]]
  allocation: Set_Inf_SndRcv() of truncation.c
  free:       Set_Inf_SndRcv() of truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Snd_MAN;

/*******************************************************
 int **Snd_GAN;
  Snd_GAN and Snd_GAN are a global atom index of which 
  informations are sent to a processor ID.
  size: Snd_GAN[numprocs][FS_Snd_Num[ID]]
  allocation: Set_Inf_SndRcv() of truncation.c
  free:       Set_Inf_SndRcv() of truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Snd_GAN;

/*******************************************************
 int **Rcv_GAN;
  Rcv_GAN are a global atom index cell index of which 
  informations are recieved at myid from a processor ID.
  size: Rcv_GAN[numprocs][F_Rcv_Num[ID]+S_Rcv_Num[ID]]
  allocation: Set_Inf_SndRcv() of truncation.c
  free:       Set_Inf_SndRcv() of truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Rcv_GAN;

/*******************************************************
 int *Num_Snd_Grid_A2B

  Num_Snd_Grid_A2B gives the number of grids data of 
  rho_i sent to ID.
  size: Num_Snd_Grid_A2B[numprocs]
  allocation: call Allocate_Arrays() in Input_std.c
  free:       call Free_Arrays in openmx.c
*******************************************************/
int *Num_Snd_Grid_A2B;

/*******************************************************
 int *Num_Rcv_Grid_A2B

  Num_Rcv_Grid_A2B gives the number of grids data of 
  rho_i received from ID.
  size: Num_Rcv_Grid_A2B[numprocs]
  allocation: call Allocate_Arrays() in Input_std.c
  free:       call Free_Arrays in openmx.c
*******************************************************/
int *Num_Rcv_Grid_A2B;

/*******************************************************
 int **Index_Snd_Grid_A2B

  Index_Snd_Grid_A2B gives indices BN, atom, and Rn 
  in the partition B associated with the grids data of 
  rho_i sent to ID.
  size: Index_Snd_Grid_A2B[numprocs][3*Num_Snd_Grid_A2B[ID]]
  allocation: allocate_grids2atoms() in truncation.c
  free:       allocate_grids2atoms() and
              call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Index_Snd_Grid_A2B;

/*******************************************************
 int **Index_Rcv_Grid_A2B

  Index_Rcv_Grid_A2B gives indices BN, atom, and Rn 
  in the partition B associated with the grids 
  data of rho_i received from ID.
  size: Index_Rcv_Grid_A2B[numprocs][3*Num_Rcv_Grid_A2B[ID]]
  allocation: allocate_grids2atoms() in truncation.c
  free:       allocate_grids2atoms() and
              call as Free_Arrays(0) in openmx.c
*******************************************************/
int **Index_Rcv_Grid_A2B;

/*******************************************************
 double ***CntCoes;
  contraction coefficients of basis orbitals
  size: CntCoes[Matomnum+MatomnumF+1]
               [List_YOUSO[7]]
               [List_YOUSO[24]]
  allocation: in truncation.c
  free:       in truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***CntCoes;

/*******************************************************
 double **atv;
  xyz translation vectors of periodically copied cells
  size: atv[TCpyCell+1][4]
  allocation: in Set_Periodic() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **atv;

/*******************************************************
 int **ratv;
  queue number of periodically copied cells
  size: ratv[TCpyCell*2+4][TCpyCell*2+4][TCpyCell*2+4];
  allocation: in Set_Periodic() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int ***ratv;

/*******************************************************
 int **atv_ijk;
  i,j,and j number of periodically copied cells
  size: atv_ijk[TCpyCell+1][4];
  allocation: in Set_Periodic() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **atv_ijk;

/*******************************************************
 double **Gxyz; 
 atomic global coordinates, velocities, and gradients of
 the total energy with respect to the atomic coordinates
  size: Gxyz[atomnum+1][YOUSO26]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **Gxyz;

/*******************************************************
 int **natn; 
  grobal index number of neighboring atoms of an atom ct_AN
  size: natn[atomnum+1][Max_FSNAN*ScaleSize+1]
  allocation: call as Allocate_Arrays(3) in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **natn;

/*******************************************************
 double ******DM;
  current and old density matrices
  size:  DM[1]
           [SpinP_switch+1]
           [Matomnum+1]
           [FNAN[Gc_AN]+1]
           [Spe_Total_NO[Cwan]]
           [Spe_Total_NO[Hwan]] 
  allocation: allocate in truncation.c
  free:       in truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
double ******DM;

/*******************************************************
 double ******iDM;
  imaginary density matrix
  size: iDM[List_YOUSO[16]]
           [2]
           [Matomnum+1]
           [FNAN[Gc_AN]+1]
           [Spe_Total_NO[Cwan]]
           [Spe_Total_NO[Hwan]] 
  allocation: allocate in truncation.c
  free:       in truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
double ******iDM;







/**************************/
/*      set in flpq       */
/**************************/

/**************************/
/* int */
/**************************/
/* allocated and determined in Set_Grid.c */
int Max_GridN_Atom,Max_OneD_Grids;
int My_NumGridB_AB;

/*******************************************************
 int ***GListTAtoms1;
  grid index (local for ct_AN) overlaping between
  two orbitals 
  size: GListTAtoms1[Matomnum+1]
                    [FNAN[Gc_AN]+1]
                    [NumOLG[Mc_AN][h_AN]]
  allocation: UCell_Box() of truncation.c
  free:       truncation.c and Free_Arrays(0) in openmx.c
*******************************************************/
int ***GListTAtoms1;

/*******************************************************
 int ***GListTAtoms2;
  grid index (local for h_AN) overlaping between
  two orbitals 
  size: GListTAtoms2[Matomnum+1]
                    [FNAN[Gc_AN]+1]
                    [NumOLG[Mc_AN][h_AN]]
  allocation: UCell_Box() of truncation.c
  free:       truncation.c and Free_Arrays(0) in openmx.c
*******************************************************/
int ***GListTAtoms2;

/*******************************************************
 int **GridListAtom; 
  neighboring grid points of an atom Mc_AN
  size: GridListAtom[Matomnum+1][List_YOUSO[11]]
  allocation: allocate in UCell_Box() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **GridListAtom;

/*******************************************************
 int **CellListAtom; 
  cell number of neighboring grid points of an atom Mc_AN
  size: CellListAtom[Matomnum+1][Max_GridN_Atom*ScaleSize+1]
  allocation: allocate in UCell_Box() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **CellListAtom;

/*******************************************************
 int **MGridListAtom; 
  neighboring grid points (medium variable) of an atom Mc_AN
  size: MGridListAtom[Matomnum+1][Max_GridN_Atom*ScaleSize+1]
  allocation: allocate in UCell_Box() of truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
int **MGridListAtom;

/*******************************************************
 int **NumOLG;
  the number of overlapping grids between atom Mc_AN 
  and atom Lh_AN
  size: NumOLG[Matomnum+1]
              [FNAN[Gc_AN]+1]
  allocation: allocate in truncation.c
  free:       in truncation.c
              and call as Free_Arrays(0) in openmx.c
*******************************************************/
int **NumOLG;

/*******************************************************
 double **Cell_Gxyz; 
 atomic global coordinates spanned
 by the unit cell vectors
  size: Cell_Gxyz[atomnum+1][4]
  allocation: call as Allocate_Arrays(1) in Input_std.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **Cell_Gxyz;


int Min_Grid_Index[4],Max_Grid_Index[4];
int Min_Grid_Index_D[4],Max_Grid_Index_D[4];

double rtv[4][4];
double gtv[4][4],rgtv[4][4],length_gtv[4];

double Grid_Origin[4];
double **Grid_Origin_N;

/**************************/
/*     type_Orbs_Grid     */
/* allocated in Set_Grid  */
/**************************/

/*******************************************************
 double **Density_Grid_B; 
  electron densities on grids in the partition B
  size: Density_Grid[2 or 4][My_NumGridB_AB]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **Density_Grid_B;

/*******************************************************
 double ***dDensity_Grid_B; 
  electron densities on grids in the partition B
  size: dDensity_Grid[3][6][My_NumGridB_AB]
       (n_d)_{ab}^k = \psi^{a*} \partial_k \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       dDensity_Grid_B[k][0][BN]　:  Re[nd_{11}^k]
       dDensity_Grid_B[k][1][BN]　:  Re[nd_{22}^k] 
       dDensity_Grid_B[k][2][BN]　:  Re[nd_{12}^k]
       dDensity_Grid_B[k][3][BN]　:  Im[nd_{12}^k]
       dDensity_Grid_B[k][4][BN]　:  Re[nd_{21}^k]
       dDensity_Grid_B[k][5][BN]　: -Im[nd_{21}^k]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***dDensity_Grid_B;


/* added by fukuda 24.Jan.2015 */
/*******************************************************
 double ***ddDensity_Grid_B; 
  electron densities on grids in the partition B
  size: ddDensity_Grid[7][6][My_NumGridB_AB]
       (n_dd)_{ab}^0 = \psi^{a*} laplacian \psi^b
       (n_dd)_{ab}^1 = \psi^{a*} \partial_x \partial_x \psi^b
       (n_dd)_{ab}^2 = \psi^{a*} \partial_y \partial_y \psi^b
       (n_dd)_{ab}^3 = \psi^{a*} \partial_z \partial_z \psi^b
       (n_dd)_{ab}^4 = \psi^{a*} \partial_x \partial_y \psi^b
       (n_dd)_{ab}^5 = \psi^{a*} \partial_y \partial_z \psi^b
       (n_dd)_{ab}^6 = \psi^{a*} \partial_z \partial_x \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       ddDensity_Grid_B[k][0][BN]　:  Re[ndd_{11}^k]
       ddDensity_Grid_B[k][1][BN]　:  Re[ndd_{22}^k] 
       ddDensity_Grid_B[k][2][BN]　:  Re[ndd_{12}^k]
       ddDensity_Grid_B[k][3][BN]　:  Im[ndd_{12}^k]
       ddDensity_Grid_B[k][4][BN]　:  Re[ndd_{21}^k]
       ddDensity_Grid_B[k][5][BN]　: -Im[ndd_{21}^k]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***ddDensity_Grid_B;

/* added by fukuda 2.May.2015 */
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
       dd2Density_Grid_B[k][0][BN]　:  Re[ndd2_{11}^k]
       dd2Density_Grid_B[k][1][BN]　:  Re[ndd2_{22}^k] 
       dd2Density_Grid_B[k][2][BN]　:  Re[ndd2_{12}^k]
       dd2Density_Grid_B[k][3][BN]　:  Im[ndd2_{12}^k]
       dd2Density_Grid_B[k][4][BN]　:  Re[ndd2_{21}^k]
       dd2Density_Grid_B[k][5][BN]　: -Im[ndd2_{21}^k]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***dd2Density_Grid_B;

/* added by fukuda 23.Jun.2015 */
/*******************************************************
 double **Density_Grid_B_i; 
  electron densities on grids in the partition B
  size: Density_Grid[2][My_NumGridB_AB]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double **Density_Grid_B_i;

/* added by fukuda 23.Jun.2015 */
/*******************************************************
 double ***dDensity_Grid_B_i; 
  electron densities on grids in the partition B
  size: dDensity_Grid[3][2][My_NumGridB_AB]
       (n_d)_{ab}^k = \psi^{a*} \partial_k \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       dDensity_Grid_B_i[k][0][BN]　:  Im[nd_{11}^k]
       dDensity_Grid_B_i[k][1][BN]　:  Im[nd_{22}^k] 
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***dDensity_Grid_B_i;


/* added by fukuda 23.Jun.2015 */
/*******************************************************
 double ***ddDensity_Grid_B_i; 
  electron densities on grids in the partition B
  size: ddDensity_Grid[7][2][My_NumGridB_AB]
       (n_dd)_{ab}^0 = \psi^{a*} laplacian \psi^b
       (n_dd)_{ab}^1 = \psi^{a*} \partial_x \partial_x \psi^b
       (n_dd)_{ab}^2 = \psi^{a*} \partial_y \partial_y \psi^b
       (n_dd)_{ab}^3 = \psi^{a*} \partial_z \partial_z \psi^b
       (n_dd)_{ab}^4 = \psi^{a*} \partial_x \partial_y \psi^b
       (n_dd)_{ab}^5 = \psi^{a*} \partial_y \partial_z \psi^b
       (n_dd)_{ab}^6 = \psi^{a*} \partial_z \partial_x \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       ddDensity_Grid_B_i[k][0][BN]　:  Im[ndd_{11}^k]
       ddDensity_Grid_B_i[k][1][BN]　:  Im[ndd_{22}^k] 
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***ddDensity_Grid_B_i;

/* added by fukuda 23.Jun.2015 */
/*******************************************************
 double ***dd2Density_Grid_B_i; 
  electron densities on grids in the partition B
  size: dd2Density_Grid[7][2][My_NumGridB_AB]
       (n_dd2)_{ab}^0 = 0
       (n_dd2)_{ab}^1 = (\partial_x \psi)^{a*} \partial_x \psi^b
       (n_dd2)_{ab}^2 = (\partial_y \psi)^{a*} \partial_y \psi^b
       (n_dd2)_{ab}^3 = (\partial_z \psi)^{a*} \partial_z \psi^b
       (n_dd2)_{ab}^4 = (\partial_x \psi)^{a*} \partial_y \psi^b
       (n_dd2)_{ab}^5 = (\partial_y \psi)^{a*} \partial_z \psi^b
       (n_dd2)_{ab}^6 = (\partial_z \psi)^{a*} \partial_x \psi^b
       -- psi^1 : upper component of psi
       -- psi^2 : lower component of psi
       dd2Density_Grid_B_i[k][0][BN]　:  Im[ndd2_{11}^k]
       dd2Density_Grid_B_i[k][1][BN]　:  Im[ndd2_{22}^k] 
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
double ***dd2Density_Grid_B_i;

/*******************************************************
 Type_Orbs_Grid ***Orbs_Grid;
  values of basis orbitals on grids
  size: Orbs_Grid[Matomnum+1]
                 [GridN_Atom[Gc_AN]]
                 [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ***Orbs_Grid;

/* added by fukuda 30.Dec.2014 */
/*******************************************************
 Type_Orbs_Grid ****dOrbs_Grid;
  values of basis orbitals on grids
  size: dOrbs_Grid[3]
                  [Matomnum+1]
                  [GridN_Atom[Gc_AN]]
                  [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ****dOrbs_Grid;

/* added by fukuda 24.Jan.2015 */
/*******************************************************
 Type_Orbs_Grid ****ddOrbs_Grid;
  values of basis orbitals on grids
  size: ddOrbs_Grid[7]
                   [Matomnum+1]
                   [GridN_Atom[Gc_AN]]
                   [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ****ddOrbs_Grid;

/*******************************************************
 Type_Orbs_Grid ***COrbs_Grid;
  values of contrated basis orbitals on grids
  size: COrbs_Grid[Matomnum+MatomnumF+1]
                  [Spe_Total_NO[Cwan]]
                  [GridN_Atom[Gc_AN]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
* ******************************************************/
Type_Orbs_Grid ***COrbs_Grid;

/* added by fukuda 30.Dec.2014 */
/*******************************************************
 Type_Orbs_Grid ****CdOrbs_Grid;
  values of contrated basis orbitals on grids
  size: CdOrbs_Grid[3]
                   [Matomnum+MatomnumF+1]
                   [Spe_Total_NO[Cwan]]
                   [GridN_Atom[Gc_AN]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ****CdOrbs_Grid;

/* added by fukuda 24.Jan.2015 */
/*******************************************************
 Type_Orbs_Grid ****CddOrbs_Grid;
  values of contrated basis orbitals on grids
  size: CddOrbs_Grid[7]
                    [Matomnum+MatomnumF+1]
                    [Spe_Total_NO[Cwan]]
                    [GridN_Atom[Gc_AN]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ****CddOrbs_Grid;

/*******************************************************
 Type_Orbs_Grid ***Orbs_Grid_FNAN;
  values of basis orbitals on grids for neighbor atoms
  which do not belong to my process.
  size: Orbs_Grid_FNAN[Matomnum+1]
                      [FNAN[Gc_AN]+1]
                      [NumOLG[Mc_AN][h_AN]]
                      [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid ****Orbs_Grid_FNAN;

/* added by fukuda 30.Dec.2014 */
/*******************************************************
 Type_Orbs_Grid ****dOrbs_Grid_FNAN;
  values of basis orbitals on grids for neighbor atoms
  which do not belong to my process.
  size: dOrbs_Grid_FNAN[3]
                       [Matomnum+1]
                       [FNAN[Gc_AN]+1]
                       [NumOLG[Mc_AN][h_AN]]
                       [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid *****dOrbs_Grid_FNAN;

/* added by fukuda 24.Jan.2015 */
/*******************************************************
 Type_Orbs_Grid ****ddOrbs_Grid_FNAN;
  values of basis orbitals on grids for neighbor atoms
  which do not belong to my process.
  size: ddOrbs_Grid_FNAN[7]
                        [Matomnum+1]
                        [FNAN[Gc_AN]+1]
                        [NumOLG[Mc_AN][h_AN]]
                        [Spe_Total_NO[Cwan]]
  allocation: allocate in truncation.c
  free:       call as Free_Arrays(0) in openmx.c
*******************************************************/
Type_Orbs_Grid *****ddOrbs_Grid_FNAN;




/********************************************/
/****     interface of subroutines       ****/
/********************************************/
void Input_std(char *file);
//void read_input_flpq();
void read_input_openmx_lpq();
void read_input_DM( int kloop, int mu );
int Set_Allocate_Atom2CPU(int weight_flag);
void Set_Inf_SndRcv();
void Set_Grid();
void Get_ddOrbitals(int wan, double x, double y, double z, double **ddChi);
void Get_Cnt_ddOrbitals(int Mc_AN, double x, double y, double z, double **ddChi);
double Set_ddOrbitals_Grid();
double Set_Density_Grid_fuku(double *****CDM);
double Set_dDensity_Grid(double *****CDM);
double Set_ddDensity_Grid(double *****CDM);
double Set_dd2Density_Grid(double *****CDM);
double Set_Density_Grid_fuku_i(double *****CDM);
double Set_dDensity_Grid_i(double *****CDM);
double Set_ddDensity_Grid_i(double *****CDM);
double Set_dd2Density_Grid_i(double *****CDM);
void out_LPQ();
void Set_Grid_Origin_N();
void Free_Arrays();
void Free_Arrays_in_Set_Grid();

void GN2N(int GN, int N3[4]);
void xyz2spherical(double x, double y, double z,
                   double xo, double yo, double zo,
                   double S_coordinate[3]);
double sgn(double nu);

void Cross_Product(double a[4], double b[4], double c[4]);
double Dot_Product(double a[4], double b[4]);
void qsort_int(long n, int *a, int *b);
void qsort_int1(long n, int *a);
void qsort_double_int(long n, double *a, int *b);
typedef struct { int a,b; } ilists;
int ilists_cmp(const ilists *x, const ilists *y);
void Get_Grid_XYZ(int GN, double xyz[4]);

void Find_CGrids(int Real_Position, int n1, int n2, int n3,
                 double Cxyz[4], int NOC[4]);
void lapack_dstegr1(INTEGER N, INTEGER EVmax, double *D, double *E, double *W, double **ev);
