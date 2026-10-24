#ifndef FLPQ_LOCAL_PHYSICAL_QUANTITIES_H
#define FLPQ_LOCAL_PHYSICAL_QUANTITIES_H

#include <stdio.h>

/* Values remain double precision until an output backend formats them.
 * Cartesian tensor entries are row-major: xx,xy,xz,yx,yy,yz,zx,zy,zz.
 * stress_diag stores three eigenvalues followed by three eigenvectors.
 */
typedef struct {
  int rank;
  int index;
  int components;
  double values[12];
} FLPQ_LocalValue;

int flpq_evaluate_local_quantity(int rank, int index, int BN,
                                 FLPQ_LocalValue *result);
void flpq_write_local_value(FILE *fp, const FLPQ_LocalValue *result,
                            int flag_cube);

void scalar_density(int n, FILE *fp, int BN, int flag_cube);
void vector_density(int n, FILE *fp, int BN);
void tensor_density(int n, FILE *fp, int BN);

#endif
