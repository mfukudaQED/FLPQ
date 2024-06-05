#include "f77func.h"

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif

#ifndef ___logical_definition___
typedef long int logical;
#define ___logical_definition___ 
#endif

#ifndef ___dcomplex_definition___
typedef struct { double r,i; } dcomplex;
#define ___dcomplex_definition___ 
#endif

void dsyevd_(char *JOBZ, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, double *W, double *WORK, 
         INTEGER *LWORK, INTEGER *IWORK,
                       INTEGER  *LIWORK, INTEGER *INFO);

void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, 
          double *VL, double *VU,
          INTEGER *IL, INTEGER *IU, double *ABSTOL, INTEGER *M, double *W, double *Z, 
           INTEGER *LDZ, double *WORK,
          INTEGER *LWORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);

int dstevx_(char *JOBZ, char *RANGE, INTEGER *N, double *D, double *E, double *VL, double *VU, INTEGER *IL, INTEGER *IU, double *ABSTOL,
                        INTEGER *M, double *W, double *Z, INTEGER *LDZ, double *WORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);

int dstegr_(char *JOBZ , char *RANGE , INTEGER *N , double *D , double *E , 
       double *VL , double *VU , INTEGER *IL , INTEGER *IU , double *ABSTOL , 
       INTEGER *M , double *W , double *Z , INTEGER *LDZ , INTEGER *ISUPPZ , 
       double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , INTEGER *INFO
       );

int dstedc_(char *COMPZ , INTEGER *N , double *D , double *E , double *Z , 
       INTEGER *LDZ , double *WORK , INTEGER *LWORK , INTEGER *IWORK , INTEGER *LIWORK , 
       INTEGER *INFO);
