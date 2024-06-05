/**********************************************************************
  util.c:

     util.c is a group of some useful subroutines.
     These subroutines are as same as ones in OpenMX.

  Log of util.c:

     01/Jul/2016  Released by M.Fukuda
***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "flpq.h"

typedef struct {
    double a,b;
} dlists;

typedef struct {
  double a;
  int b;
} dilists;

static void Find_ln(int i0, int n0, int res[2]);
static int dlists_cmp(const dlists *x, const dlists *y);


void GN2N(int GN, int N3[4])
{
  int n1,n2,n3;

  n1 = GN/(Ngrid2*Ngrid3);
  n2 = (GN - n1*(Ngrid2*Ngrid3))/Ngrid3;
  n3 = GN - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;
  N3[1] = n1;
  N3[2] = n2;
  N3[3] = n3;
}


void xyz2spherical(double x, double y, double z,
                   double xo, double yo, double zo,
                   double S_coordinate[3])
{
  double dx,dy,dz,r,r1,theta,phi,dum,dum1,Min_r;

  Min_r = 10e-15;

  dx = x - xo;
  dy = y - yo;
  dz = z - zo;

  dum = dx*dx + dy*dy; 
  r = sqrt(dum + dz*dz);
  r1 = sqrt(dum);

  if (Min_r<=r){

    if (r<fabs(dz))
      dum1 = sgn(dz)*1.0;
    else
      dum1 = dz/r;

    theta = acos(dum1);

    if (Min_r<=r1){
      if (0.0<=dx){

        if (r1<fabs(dy))
          dum1 = sgn(dy)*1.0;
        else
          dum1 = dy/r1;        
  
        phi = asin(dum1);
      }
      else{

        if (r1<fabs(dy))
          dum1 = sgn(dy)*1.0;
        else
          dum1 = dy/r1;        

        phi = PI - asin(dum1);
      }
    }
    else{
      phi = 0.0;
    }
  }
  else{
    theta = 0.5*PI;
    phi = 0.0;
  }

  S_coordinate[0] = r;
  S_coordinate[1] = theta;
  S_coordinate[2] = phi;
}

double sgn(double nu)
{
  double result;
  if (nu<0.0)
    result = -1.0;
  else
    result = 1.0;
  return result;
}


void Cross_Product(double a[4], double b[4], double c[4])
{
  c[1] = a[2]*b[3] - a[3]*b[2]; 
  c[2] = a[3]*b[1] - a[1]*b[3]; 
  c[3] = a[1]*b[2] - a[2]*b[1];
}

double Dot_Product(double a[4], double b[4])
{
  double sum;
  sum = a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; 
  return sum;
}


void qsort_int(long n, int *a, int *b)
{
  int i;
  ilists *AB;

  AB = (ilists*)malloc(sizeof(ilists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i+1];
    AB[i].b = b[i+1];
  }

  qsort(AB, n, sizeof(ilists), (int(*)(const void*, const void*))ilists_cmp);

  for (i=0; i<n; i++){
    a[i+1] = AB[i].a;
    b[i+1] = AB[i].b;
  }

  free(AB);
}

void qsort_int1(long n, int *a)
{
  qsort(a, n, sizeof(int), (int(*)(const void*, const void*))ilists_cmp);
}

void qsort_double_int(long n, double *a, int *b)
{
  int i;
  dilists *AB;

  AB = (dilists*)malloc(sizeof(dilists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i];     
    AB[i].b = b[i];
  }

  qsort(AB, n, sizeof(dilists), (int(*)(const void*, const void*))dlists_cmp);

  for (i=0; i<n; i++){
    a[i] = AB[i].a;
    b[i] = AB[i].b;
  }

  free(AB);
}
 
int dlists_cmp(const dlists *x, const dlists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}

int ilists_cmp(const ilists *x, const ilists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}


void Get_Grid_XYZ(int GN, double xyz[4])
{
  int n1,n2,n3;

  n1 = GN/(Ngrid2*Ngrid3);
  n2 = (GN - n1*(Ngrid2*Ngrid3))/Ngrid3;
  n3 = GN - n1*(Ngrid2*Ngrid3) - n2*Ngrid3;

  xyz[1] = (double)n1*gtv[1][1] + (double)n2*gtv[2][1]
         + (double)n3*gtv[3][1] + Grid_Origin[1];
  xyz[2] = (double)n1*gtv[1][2] + (double)n2*gtv[2][2]
         + (double)n3*gtv[3][2] + Grid_Origin[2];
  xyz[3] = (double)n1*gtv[1][3] + (double)n2*gtv[2][3]
         + (double)n3*gtv[3][3] + Grid_Origin[3];
}


void Find_CGrids(int Real_Position, int n1, int n2, int n3,
                 double Cxyz[4], int NOC[4])
{
  /*********************************************************
   Real_Position==0
       gives the coordinates in the original cell (Rn==0)

   Real_Position==1
       gives the coordinates in the translated cell
  **********************************************************/

  int Rn,l1,l2,l3,N;
  int nn1,nn2,nn3,res[2];
  double x0,y0,z0;
  
  Find_ln(1,n1,res);
  l1  = res[0];
  nn1 = res[1];
  
  Find_ln(2,n2,res);
  l2  = res[0];
  nn2 = res[1];
  
  Find_ln(3,n3,res);
  l3  = res[0];
  nn3 = res[1];
   
  if (CpyCell<abs(l1) || CpyCell<abs(l2) || CpyCell<abs(l3)){

    /* outside of Tcell */

    Rn = 0;
    x0 = atv[Rn][1];
    y0 = atv[Rn][2];
    z0 = atv[Rn][3];

    Cxyz[1] = 10e+5;
    Cxyz[2] = 10e+5;
    Cxyz[3] = 10e+5;

  }
  else{
  
    //Rn = R_atv(CpyCell,l1,l2,l3);
    Rn = ratv[l1+CpyCell][l2+CpyCell][l3+CpyCell];
    x0 = atv[Rn][1];
    y0 = atv[Rn][2];
    z0 = atv[Rn][3];
  
    /*  N = nn1*Ngrid2*Ngrid3 + nn2*Ngrid3 + nn3;  */
  
    if (Real_Position==0){
      Cxyz[1] = (double)nn1*gtv[1][1] + (double)nn2*gtv[2][1]
              + (double)nn3*gtv[3][1] + Grid_Origin[1];
      Cxyz[2] = (double)nn1*gtv[1][2] + (double)nn2*gtv[2][2]
              + (double)nn3*gtv[3][2] + Grid_Origin[2];
      Cxyz[3] = (double)nn1*gtv[1][3] + (double)nn2*gtv[2][3]
              + (double)nn3*gtv[3][3] + Grid_Origin[3];
    }
    else{
      Cxyz[1] = x0 
              + (double)nn1*gtv[1][1] + (double)nn2*gtv[2][1]
              + (double)nn3*gtv[3][1] + Grid_Origin[1];
      Cxyz[2] = y0
              + (double)nn1*gtv[1][2] + (double)nn2*gtv[2][2]
              + (double)nn3*gtv[3][2] + Grid_Origin[2];
      Cxyz[3] = z0
              + (double)nn1*gtv[1][3] + (double)nn2*gtv[2][3]
              + (double)nn3*gtv[3][3] + Grid_Origin[3];
    }
  }

  NOC[0] = Rn;
  NOC[1] = nn1;
  NOC[2] = nn2;
  NOC[3] = nn3;
}


void Find_ln(int i0, int n0, int res[2])
{
  int l0,n1,N;

  if      (i0==1) N = Ngrid1;
  else if (i0==2) N = Ngrid2;
  else if (i0==3) N = Ngrid3;

  if (n0<0){
    l0 = -((abs(n0)-1)/N + 1);
    n1 = N - (abs(n0) - N*(abs(l0)-1));
  }  
  else if (N<=n0){
    l0 = n0/N;
    n1 = n0 - N*l0;
  }
  else{
    l0 = 0;  
    n1 = n0;
  }

  res[0] = l0;
  res[1] = n1; 
}

