/**********************************************************************
  Get_Cnt_ddOrbitals.c:

     Get_Cnt_ddOrbitals.c is a subrutine to calculate
     derivatives of contracted basis orbitals

  Log of Get_Cnt_ddOrbitals.c:

     23/Jan/2015  Released by M.Fukuda

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "flpq.h"

void Get_Cnt_ddOrbitals(int Mc_AN, double x, double y, double z, double **ddChi)
{
  static int firsttime=1;
  int Gc_AN,i,j,L0,Mul0,M0,i1,al,wan,p;
  int mp_min,mp_max,m,po;
  double S_coordinate[3];
  double dum,dum1,dum2,dum3,dum4;
  double siQ,coQ,siP,coP,Q,P,R,Rmin;
  double **RF;
  double **dRF;
  double **ddRF;
  double **lapRF;
  double **AF;
  double **dAFQ;
  double **dAFP;
  double **ddAFQQ;
  double **ddAFQP;
  double **ddAFPP;
  double h1,h2,h3,f1,f2,f3,f4,a,b;
  double g1,g2,x1,x2,y1,y2,y12,y22,f,df;
  double ddf;
  double dRx,dRy,dRz,dQx,dQy,dQz,dPx,dPy,dPz;
  double ddRxx,ddRyy,ddRzz,ddRxy,ddRyz,ddRzx;
  double ddQxx,ddQyy,ddQzz,ddQxy,ddQyz,ddQzx;
  double ddPxx,ddPyy,ddPzz,ddPxy,ddPyz,ddPzx;
  double dChiR,dChiQ,dChiP,tmp0,tmp1,tmp2,tmp3;
  double SH[Supported_MaxL*2+1][2];
  double dSHt[Supported_MaxL*2+1][2];
  double dSHp[Supported_MaxL*2+1][2];
  double dSHtt[Supported_MaxL*2+1][2];
  double dSHtp[Supported_MaxL*2+1][2];
  double dSHpp[Supported_MaxL*2+1][2];
  double Rw,Rw2;
  double dRxdRx,dRydRy,dRzdRz,dRxdRy,dRydRz,dRzdRx;
  double dPxdPx,dPydPy,dPzdPz,dPxdPy,dPydPz,dPzdPx;
  double dQxdQx,dQydQy,dQzdQz,dQxdQy,dQydQz,dQzdQx;
  double dQxdPx,dQydPy,dQzdPz,dQxdPy,dQydPz,dQzdPx;
  double dQydPx,dQzdPy,dQxdPz;
  double dQxdPy_dQydPx, dQydPz_dQzdPy, dQzdPx_dQxdPz;
  double dRFx,dRFy,dRFz,ddRFxx,ddRFyy,ddRFzz,ddRFxy,ddRFyz,ddRFzx;
  double dAFx,dAFy,dAFz,ddAFxx,ddAFyy,ddAFzz,ddAFxy,ddAFyz,ddAFzx;

  /****************************************************
   allocation of arrays:
  
   double  RF[List_YOUSO[25]+1][List_YOUSO[24]];
   double dRF[List_YOUSO[25]+1][List_YOUSO[24]];
   double   AF[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
   double dAFQ[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
   double dAFP[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
  ****************************************************/
  
  RF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    RF[i] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
  }

  dRF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    dRF[i] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
  }

  ddRF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    ddRF[i] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
  }

  lapRF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    lapRF[i] = (double*)malloc(sizeof(double)*List_YOUSO[24]);
  }

  AF = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    AF[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  dAFQ = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    dAFQ[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  dAFP = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    dAFP[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  ddAFQQ = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    ddAFQQ[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  ddAFQP = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    ddAFQP[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  ddAFPP = (double**)malloc(sizeof(double*)*(List_YOUSO[25]+1));
  for (i=0; i<(List_YOUSO[25]+1); i++){
    ddAFPP[i] = (double*)malloc(sizeof(double)*(2*(List_YOUSO[25]+1)+1));
  }

  /* start calc. */

  Gc_AN = M2G[Mc_AN];
  //Rmin = 10e-12;
  Rmin = 10e-14;
  wan = WhatSpecies[Gc_AN];
  xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate);
  R = S_coordinate[0];
  Q = S_coordinate[1];
  P = S_coordinate[2];

  if (R<Rmin){
    x = x + Rmin;
    y = y + Rmin;
    z = z + Rmin;
    xyz2spherical(x,y,z,0.0,0.0,0.0,S_coordinate);
    R = S_coordinate[0];
    Q = S_coordinate[1];
    P = S_coordinate[2];
  }  

  po = 0;
  mp_min = 0;
  mp_max = Spe_Num_Mesh_PAO[wan] - 1;

  if (Spe_PAO_RV[wan][Spe_Num_Mesh_PAO[wan]-1]<R){
    for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){
         RF[L0][Mul0] = 0.0;
        dRF[L0][Mul0] = 0.0;
       ddRF[L0][Mul0] = 0.0;
       lapRF[L0][Mul0] = 0.0;
      }
    }
    po = 1;
  }
  else if (R<Spe_PAO_RV[wan][0]){

    h1 = Spe_PAO_RV[wan][0] - Spe_PAO_RV[wan][1];

    for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){

        y1 = Spe_PAO_RWF[wan][L0][Mul0][0] - Spe_PAO_RWF[wan][L0][Mul0][1];
        a = y1/h1;
        b = Spe_PAO_RWF[wan][L0][Mul0][0] - a*Spe_PAO_RV[wan][0];

         RF[L0][Mul0] = a*R + b;
        dRF[L0][Mul0] = a;
       ddRF[L0][Mul0] = 0.0;
       lapRF[L0][Mul0] = 0.0;
      }
    }

  }
  else{

    do{
      m = (mp_min + mp_max)/2;
      if (Spe_PAO_RV[wan][m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    h1 = Spe_PAO_RV[wan][m-1] - Spe_PAO_RV[wan][m-2];
    h2 = Spe_PAO_RV[wan][m]   - Spe_PAO_RV[wan][m-1];
    h3 = Spe_PAO_RV[wan][m+1] - Spe_PAO_RV[wan][m];

    x1 = R - Spe_PAO_RV[wan][m-1];
    x2 = R - Spe_PAO_RV[wan][m];
    y1 = x1/h2;
    y2 = x2/h2;
    y12 = y1*y1;
    y22 = y2*y2;

    dum = h1 + h2;
    dum1 = h1/h2/dum;
    dum2 = h2/h1/dum;
    dum = h2 + h3;
    dum3 = h2/h3/dum;
    dum4 = h3/h2/dum;

    for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
      for (Mul0=0; Mul0<Spe_Num_Basis[wan][L0]; Mul0++){

        f1 = Spe_PAO_RWF[wan][L0][Mul0][m-2];
        f2 = Spe_PAO_RWF[wan][L0][Mul0][m-1];
        f3 = Spe_PAO_RWF[wan][L0][Mul0][m];
        f4 = Spe_PAO_RWF[wan][L0][Mul0][m+1];

        if (m==1){
          h1 = -(h2+h3);
          f1 = f4;
        }
        else if (m==(Spe_Num_Mesh_PAO[wan]-1)){
          h3 = -(h1+h2);
          f4 = f1;
        }

        dum = f3 - f2;
        g1 = dum*dum1 + (f2-f1)*dum2;
        g2 = (f4-f3)*dum3 + dum*dum4;

        f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
           + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

        df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
           + y22*(2.0*f2 + h2*g1)/h2
           + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
           - y12*(2.0*f3 - h2*g2)/h2;
        
        ddf = 2.0/h2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
            + 4.0*y2*(2.0*f2 + h2*g1)/h2/h2
            + 2.0/h2/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
            - 4.0*y1*(2.0*f3 - h2*g2)/h2/h2;

          RF[L0][Mul0] = f;
         dRF[L0][Mul0] = df;
        ddRF[L0][Mul0] = ddf;
       lapRF[L0][Mul0] = (2.0*R*dRF[L0][Mul0] + R*R*ddRF[L0][Mul0] - (double)(L0*(L0+1))*RF[L0][Mul0])/(R*R);

      }
    } 
  }

  /* dr/dx,y,z, dQ/dx,y,z, dP/dx,y,z and dAngular */

  if (po==0){

    /* Angular */
    siQ = sin(Q);
    coQ = cos(Q);
    siP = sin(P);
    coP = cos(P);
    Rw = 1.0/R;
    Rw2 = Rw*Rw;

    dRx = siQ*coP;
    dRy = siQ*siP;
    dRz = coQ;

    ddRxx = (1.0 - siQ*siQ*coP*coP)*Rw;
    ddRyy = (1.0 - siQ*siQ*siP*siP)*Rw;
    ddRzz = (1.0 - coQ*coQ)*Rw;
    ddRxy = -(siQ*siQ*coP*siP)*Rw;
    ddRyz = -(coQ*siQ*siP)*Rw;
    ddRzx = -(coQ*siQ*coP)*Rw;

    if (Rmin<R){
      dQx = coQ*coP*Rw;
      dQy = coQ*siP*Rw;
      dQz = -siQ*Rw;

      ddQxx = (-2.0*coQ*siQ*coP*coP + (coQ*siP*siP)/siQ)*Rw2;
      ddQyy = (-2.0*coQ*siQ*siP*siP + (coQ*coP*coP)/siQ)*Rw2;
      ddQzz = 2.0*coQ*siQ*Rw2;
      ddQxy = -(2.0*coQ*siQ*coP*siP + (coQ*coP*siP)/siQ)*Rw2;
      ddQyz = -(coQ*coQ - siQ*siQ)*siP*Rw2;
      ddQzx = -(coQ*coQ - siQ*siQ)*coP*Rw2;
    }
    else{
      dQx = 0.0;
      dQy = 0.0;
      dQz = 0.0;
      ddQxx = 0.0; 
      ddQyy = 0.0;
      ddQzz = 0.0;
      ddQxy = 0.0;
      ddQyz = 0.0;
      ddQzx = 0.0;
    }

    /* RICS note 72P */

    if (Rmin<R){
      dPx = -siP/siQ*Rw;
      dPy = coP/siQ*Rw;
      dPz = 0.0;
      ddPxx =  2.0*coP*siP/siQ/siQ*Rw2;
      ddPyy = -2.0*coP*siP/siQ/siQ*Rw2;
      ddPzz = 0.0;
      ddPxy = -(coP*coP - siP*siP)/siQ/siQ*Rw2;
      ddPyz = 0.0;
      ddPzx = 0.0;
    }
    else{
      dPx = 0.0;
      dPy = 0.0;
      dPz = 0.0;
      ddPxx = 0.0;
      ddPyy = 0.0;
      ddPzz = 0.0;
      ddPxy = 0.0;
      ddPyz = 0.0;
      ddPzx = 0.0;
    }

    for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
      if (L0==0){
        AF[0][0] = 0.282094791773878;
        dAFQ[0][0] = 0.0;
        dAFP[0][0] = 0.0;
        ddAFQQ[0][0] = 0.0;
        ddAFQP[0][0] = 0.0;
        ddAFPP[0][0] = 0.0;
      }
      else if (L0==1){
        dum = 0.48860251190292*siQ;

        AF[1][0] = dum*coP;
        AF[1][1] = dum*siP;
        AF[1][2] = 0.48860251190292*coQ;

        dAFQ[1][0] = 0.48860251190292*coQ*coP;
        dAFQ[1][1] = 0.48860251190292*coQ*siP;
        dAFQ[1][2] = -0.48860251190292*siQ;

        dAFP[1][0] = -0.48860251190292*siQ*siP;
        dAFP[1][1] = 0.48860251190292*siQ*coP;
        dAFP[1][2] = 0.0;

        ddAFQQ[1][0] = -0.48860251190292*siQ*coP;
        ddAFQQ[1][1] = -0.48860251190292*siQ*siP;
        ddAFQQ[1][2] = -0.48860251190292*coQ;

        ddAFPP[1][0] = -0.48860251190292*siQ*coP;
        ddAFPP[1][1] = -0.48860251190292*siQ*siP;
        ddAFPP[1][2] = 0.0;

        ddAFQP[1][0] = -0.48860251190292*coQ*siP;
        ddAFQP[1][1] = 0.48860251190292*coQ*coP;
        ddAFQP[1][2] = 0.0;
      }
      else if (L0==2){
        dum1 = siQ*siQ;
        dum2 = 1.09254843059208*siQ*coQ;
        AF[2][0] = 0.94617469575756*coQ*coQ - 0.31539156525252;
        AF[2][1] = 0.54627421529604*dum1*(1.0 - 2.0*siP*siP);
        AF[2][2] = 1.09254843059208*dum1*siP*coP;
        AF[2][3] = dum2*coP;
        AF[2][4] = dum2*siP;

        dAFQ[2][0] = -1.89234939151512*siQ*coQ;
        dAFQ[2][1] = 1.09254843059208*siQ*coQ*(1.0 - 2.0*siP*siP);
        dAFQ[2][2] = 2.18509686118416*siQ*coQ*siP*coP;
        dAFQ[2][3] = 1.09254843059208*(1.0 - 2.0*siQ*siQ)*coP;
        dAFQ[2][4] = 1.09254843059208*(1.0 - 2.0*siQ*siQ)*siP;

        /* RICS note 72P */

        dAFP[2][0] = 0.0;
        dAFP[2][1] = -2.18509686118416*siQ*siQ*siP*coP;
        dAFP[2][2] = 1.09254843059208*siQ*siQ*(1.0 - 2.0*siP*siP);
        dAFP[2][3] = -1.09254843059208*siQ*coQ*siP;
        dAFP[2][4] = 1.09254843059208*siQ*coQ*coP;

        ddAFQQ[2][0] = -1.89234939151512*(coQ*coQ-siQ*siQ);
        ddAFQQ[2][1] = 1.09254843059208*(coQ*coQ-siQ*siQ)*(1.0 - 2.0*siP*siP);
        ddAFQQ[2][2] = 2.18509686118416*(coQ*coQ-siQ*siQ)*siP*coP;
        ddAFQQ[2][3] = -1.09254843059208*4.0*siQ*coQ*coP;
        ddAFQQ[2][4] = -1.09254843059208*4.0*siQ*coQ*siP;

        ddAFPP[2][0] = 0.0;
        ddAFPP[2][1] = -2.18509686118416*siQ*siQ*(coP*coP-siP*siP);
        ddAFPP[2][2] = -1.09254843059208*siQ*siQ*4.0*siP*coP;
        ddAFPP[2][3] = -1.09254843059208*siQ*coQ*coP;
        ddAFPP[2][4] = -1.09254843059208*siQ*coQ*siP;

        ddAFQP[2][0] = 0.0;
        ddAFQP[2][1] = -1.09254843059208*siQ*coQ*4.0*siP*coP;
        ddAFQP[2][2] = 2.18509686118416*siQ*coQ*(coP*coP-siP*siP);
        ddAFQP[2][3] = -1.09254843059208*(1.0 - 2.0*siQ*siQ)*siP;
        ddAFQP[2][4] = 1.09254843059208*(1.0 - 2.0*siQ*siQ)*coP;
      }
      else if (L0==3){
        AF[3][0] = 0.373176332590116*(5.0*coQ*coQ*coQ - 3.0*coQ);
        AF[3][1] = 0.457045799464466*coP*siQ*(5.0*coQ*coQ - 1.0);
        AF[3][2] = 0.457045799464466*siP*siQ*(5.0*coQ*coQ - 1.0);
        AF[3][3] = 1.44530572132028*siQ*siQ*coQ*(coP*coP-siP*siP);
        AF[3][4] = 2.89061144264055*siQ*siQ*coQ*siP*coP;
        AF[3][5] = 0.590043589926644*siQ*siQ*siQ*(4.0*coP*coP*coP - 3.0*coP);
        AF[3][6] = 0.590043589926644*siQ*siQ*siQ*(3.0*siP - 4.0*siP*siP*siP);

        dAFQ[3][0] = 0.373176332590116*siQ*(-15.0*coQ*coQ + 3.0);
        dAFQ[3][1] = 0.457045799464466*coP*coQ*(15.0*coQ*coQ - 11.0);
        dAFQ[3][2] = 0.457045799464466*siP*coQ*(15.0*coQ*coQ - 11.0);
        dAFQ[3][3] = 1.44530572132028*(coP*coP-siP*siP)*siQ*(2.0*coQ*coQ-siQ*siQ);
        dAFQ[3][4] = 2.89061144264055*coP*siP*siQ*(2.0*coQ*coQ - siQ*siQ);
        dAFQ[3][5] = 1.770130769779932*coP*coQ*siQ*siQ*(-3.0 + 4.0*coP*coP);
        dAFQ[3][6] = 1.770130769779932*coQ*siP*siQ*siQ*( 3.0 - 4.0*siP*siP);

        /* RICS note 72P */

        dAFP[3][0] = 0.0;
        dAFP[3][1] = 0.457045799464466*siP*siQ*(-5.0*coQ*coQ + 1.0);
        dAFP[3][2] = 0.457045799464466*coP*siQ*( 5.0*coQ*coQ - 1.0);
        dAFP[3][3] = -5.781222885281120*coP*coQ*siP*siQ*siQ;
        dAFP[3][4] = 2.89061144264055*siQ*coQ*siQ*(coP*coP - siP*siP);
        dAFP[3][5] = 1.770130769779932*siP*siQ*siQ*siQ*(1.0 - 4.0*coP*coP);
        dAFP[3][6] = 1.770130769779932*coP*siQ*siQ*siQ*(1.0 - 4.0*siP*siP);

        ddAFQQ[3][0] = 0.373176332590116*3.0*coQ*(-15.0*coQ*coQ + 11.0);
        ddAFQQ[3][1] = 0.457045799464466*coP*siQ*(-45.0*coQ*coQ + 11.0);
        ddAFQQ[3][2] = 0.457045799464466*siP*siQ*(-45.0*coQ*coQ + 11.0);
        ddAFQQ[3][3] = 1.44530572132028*(coP*coP-siP*siP)*coQ*(9.0*coQ*coQ-7.0);
        ddAFQQ[3][4] = 2.89061144264055*coP*siP          *coQ*(9.0*coQ*coQ-7.0);
        ddAFQQ[3][5] = 1.770130769779932*coP*(-3.0 + 4.0*coP*coP)*siQ*(2.0*coQ*coQ-siQ*siQ);
        ddAFQQ[3][6] = 1.770130769779932*siP*( 3.0 - 4.0*siP*siP)*siQ*(2.0*coQ*coQ-siQ*siQ);

        ddAFPP[3][0] = 0.0;
        ddAFPP[3][1] =  0.457045799464466*coP*siQ*(-5.0*coQ*coQ + 1.0);
        ddAFPP[3][2] = -0.457045799464466*siP*siQ*( 5.0*coQ*coQ - 1.0);
        ddAFPP[3][3] = -5.781222885281120*(coP*coP-siP*siP)*coQ*siQ*siQ;
        ddAFPP[3][4] = -5.781222885281120*siQ*coQ*siQ*2.0*coP*siP;
        ddAFPP[3][5] = 1.770130769779932*siQ*siQ*siQ*3.0*coP*( 4.0*siP*siP-1.0);
        ddAFPP[3][6] = 1.770130769779932*siQ*siQ*siQ*3.0*siP*(-4.0*coP*coP+1.0);

        ddAFQP[3][0] = 0.0;
        ddAFQP[3][1] = -0.457045799464466*siP*coQ*(15.0*coQ*coQ - 11.0);
        ddAFQP[3][2] = 0.457045799464466*coP*coQ*(15.0*coQ*coQ - 11.0);
        ddAFQP[3][3] = -5.781222885281120*coP*siP*siQ*(2.0*coQ*coQ-siQ*siQ);
        ddAFQP[3][4] = 2.89061144264055*(coP*coP-siP*siP)*siQ*(2.0*coQ*coQ - siQ*siQ);
        ddAFQP[3][5] = 1.770130769779932*3.0*siQ*siQ*coQ*siP*(1.0 - 4.0*coP*coP);
        ddAFQP[3][6] = 1.770130769779932*3.0*siQ*siQ*coQ*coP*(1.0 - 4.0*siP*siP);
      }

      else if (4<=L0){
        printf("The g orbital is not available\n");
      }
      //else if (4<=L0){

      //  for(m=-L0; m<=L0; m++){ 
      //    //ComplexSH(L0,m,Q,P,SH[L0+m],dSHt[L0+m],dSHp[L0+m]);
      //    ComplexSH2(L0,m,Q,P,SH[L0+m],dSHt[L0+m],dSHp[L0+m],dSHtt[L0+m],dSHpp[L0+m],dSHtp[L0+m]);
      //  }

      //  AF[L0][0]   =   -SH[L0][0];
      //  dAFQ[L0][0] = -dSHt[L0][0];
      //  dAFP[L0][0] = -dSHp[L0][0];
      //  ddAFQQ[L0][0] = -dSHtt[L0][0];
      //  ddAFQP[L0][0] = -dSHtp[L0][0];
      //  ddAFPP[L0][0] = -dSHpp[L0][0];

      //  j = -1;
      //  for (i=1; i<(L0*2+1); i=i+4){
      //    j++;
      //    AF[L0][i] = -2.0/sqrt(2.0)*SH[L0-(2*j+1)][0];
      //    dAFQ[L0][i] = -2.0/sqrt(2.0)*dSHt[L0-(2*j+1)][0];
      //    dAFP[L0][i] = -2.0/sqrt(2.0)*dSHp[L0-(2*j+1)][0];
      //    ddAFQQ[L0][i] = -2.0/sqrt(2.0)*dSHtt[L0-(2*j+1)][0];
      //    ddAFQP[L0][i] = -2.0/sqrt(2.0)*dSHtp[L0-(2*j+1)][0];
      //    ddAFPP[L0][i] = -2.0/sqrt(2.0)*dSHpp[L0-(2*j+1)][0];
      //  }

      //  j = 0;
      //  for (i=3; i<(L0*2+1); i=i+4){
      //    j++;
      //    AF[L0][i] = -2.0/sqrt(2.0)*SH[L0-2*j][0];
      //    dAFQ[L0][i] = -2.0/sqrt(2.0)*dSHt[L0-2*j][0];
      //    dAFP[L0][i] = -2.0/sqrt(2.0)*dSHp[L0-2*j][0];
      //    ddAFQQ[L0][i] = -2.0/sqrt(2.0)*dSHtt[L0-2*j][0];
      //    ddAFQP[L0][i] = -2.0/sqrt(2.0)*dSHtp[L0-2*j][0];
      //    ddAFPP[L0][i] = -2.0/sqrt(2.0)*dSHpp[L0-2*j][0];
      //  }

      //  j = -1;
      //  for (i=2; i<(L0*2+1); i=i+4){
      //    j++;
      //    AF[L0][i] = -2.0/sqrt(2.0)*SH[L0-(2*j+1)][1];
      //    dAFQ[L0][i] = -2.0/sqrt(2.0)*dSHt[L0-(2*j+1)][1];
      //    dAFP[L0][i] = -2.0/sqrt(2.0)*dSHp[L0-(2*j+1)][1];
      //    ddAFQQ[L0][i] = -2.0/sqrt(2.0)*dSHtt[L0-(2*j+1)][1];
      //    ddAFQP[L0][i] = -2.0/sqrt(2.0)*dSHtp[L0-(2*j+1)][1];
      //    ddAFPP[L0][i] = -2.0/sqrt(2.0)*dSHpp[L0-(2*j+1)][1];
      //  }

      //  j = 0;
      //  for (i=4; i<(L0*2+1); i=i+4){
      //    j++;
      //    AF[L0][i] = -2.0/sqrt(2.0)*SH[L0-2*j][1];
      //    dAFQ[L0][i] = -2.0/sqrt(2.0)*dSHt[L0-2*j][1];
      //    dAFP[L0][i] = -2.0/sqrt(2.0)*dSHp[L0-2*j][1];
      //    ddAFQQ[L0][i] = -2.0/sqrt(2.0)*dSHtt[L0-2*j][1];
      //    ddAFQP[L0][i] = -2.0/sqrt(2.0)*dSHtp[L0-2*j][1];
      //    ddAFPP[L0][i] = -2.0/sqrt(2.0)*dSHpp[L0-2*j][1];
      //  }
      //}

    }
  }

  /* Contracted Chi */
  al = -1;
  for (L0=0; L0<=Spe_MaxL_Basis[wan]; L0++){
    for (Mul0=0; Mul0<Spe_Num_CBasis[wan][L0]; Mul0++){
      for (M0=0; M0<=2*L0; M0++){
        al++;

        tmp0 = 0.0;
        tmp1 = 0.0;
        tmp2 = 0.0;
        for (p=0; p<Spe_Specified_Num[wan][al]; p++){
          tmp0 = tmp0 + CntCoes[Mc_AN][al][p]*RF[L0][p];
          tmp1 = tmp1 + CntCoes[Mc_AN][al][p]*dRF[L0][p];
          tmp2 = tmp2 + CntCoes[Mc_AN][al][p]*ddRF[L0][p];
          tmp3 = tmp3 + CntCoes[Mc_AN][al][p]*lapRF[L0][p];
        }        

        dChiR = tmp1*AF[L0][M0];
        dChiQ = tmp0*dAFQ[L0][M0];
        dChiP = tmp0*dAFP[L0][M0];

        dRFx = dRx*tmp1;
        dRFy = dRy*tmp1;
        dRFz = dRz*tmp1;
        ddRFxx = ddRxx*tmp1 + dRxdRx*tmp2;
        ddRFyy = ddRyy*tmp1 + dRydRy*tmp2;
        ddRFzz = ddRzz*tmp1 + dRzdRz*tmp2;
        ddRFxy = ddRxy*tmp1 + dRxdRy*tmp2;
        ddRFyz = ddRyz*tmp1 + dRydRz*tmp2;
        ddRFzx = ddRzx*tmp1 + dRzdRx*tmp2;

        dAFx = dQx*dAFQ[L0][M0] + dPx*dAFP[L0][M0];
        dAFy = dQy*dAFQ[L0][M0] + dPy*dAFP[L0][M0];
        dAFz = dQz*dAFQ[L0][M0] + dPz*dAFP[L0][M0];
        ddAFxx = ddQxx*dAFQ[L0][M0] + ddPxx*dAFP[L0][M0] 
                + dQxdQx*ddAFQQ[L0][M0] + 2.0*dQxdPx*ddAFQP[L0][M0] + dPxdPx*ddAFPP[L0][M0];
        ddAFyy = ddQyy*dAFQ[L0][M0] + ddPyy*dAFP[L0][M0] 
                + dQydQy*ddAFQQ[L0][M0] + 2.0*dQydPy*ddAFQP[L0][M0] + dPydPy*ddAFPP[L0][M0];
        ddAFzz = ddQzz*dAFQ[L0][M0] + ddPzz*dAFP[L0][M0] 
                + dQzdQz*ddAFQQ[L0][M0] + 2.0*dQzdPz*ddAFQP[L0][M0] + dPzdPz*ddAFPP[L0][M0];
        ddAFxy = ddQxy*dAFQ[L0][M0] + ddPxy*dAFP[L0][M0] 
                + dQxdQy*ddAFQQ[L0][M0] + dQxdPy_dQydPx*ddAFQP[L0][M0] + dPxdPy*ddAFPP[L0][M0];
        ddAFyz = ddQyz*dAFQ[L0][M0] + ddPyz*dAFP[L0][M0] 
                + dQydQz*ddAFQQ[L0][M0] + dQydPz_dQzdPy*ddAFQP[L0][M0] + dPydPz*ddAFPP[L0][M0];
        ddAFzx = ddQzx*dAFQ[L0][M0] + ddPzx*dAFP[L0][M0] 
                + dQzdQx*ddAFQQ[L0][M0] + dQzdPx_dQxdPz*ddAFQP[L0][M0] + dPzdPx*ddAFPP[L0][M0];
 
        ddChi[0][al] = tmp0*AF[L0][M0];
        ddChi[1][al] = dRx*dChiR + dQx*dChiQ + dPx*dChiP;
        ddChi[2][al] = dRy*dChiR + dQy*dChiQ + dPy*dChiP;
        ddChi[3][al] = dRz*dChiR + dQz*dChiQ + dPz*dChiP;
        ddChi[4][al] = tmp3*AF[L0][M0];
        ddChi[5][i1] = ddRFxx*AF[L0][M0] + 2.0*dRFx*dAFx + tmp0*ddAFxx; /* d2/dx2 */
        ddChi[6][i1] = ddRFyy*AF[L0][M0] + 2.0*dRFy*dAFy + tmp0*ddAFyy; /* d2/dy2 */
        ddChi[7][i1] = ddRFzz*AF[L0][M0] + 2.0*dRFz*dAFz + tmp0*ddAFzz; /* d2/dz2 */
        ddChi[8][i1] = ddRFxy*AF[L0][M0] + dRFy*dAFx + dRFx*dAFy + tmp0*ddAFxy; /* d2/dxdy */
        ddChi[9][i1] = ddRFyz*AF[L0][M0] + dRFz*dAFy + dRFy*dAFz + tmp0*ddAFyz; /* d2/dydz */
        ddChi[10][i1]= ddRFzx*AF[L0][M0] + dRFx*dAFz + dRFz*dAFx + tmp0*ddAFzx; /* d2/dzdx */
      }
    }
  }

  /****************************************************
   freeing of arrays:

   double  RF[List_YOUSO[25]+1][List_YOUSO[24]];
   double dRF[List_YOUSO[25]+1][List_YOUSO[24]];
   double   AF[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
   double dAFQ[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
   double dAFP[List_YOUSO[25]+1][2*(List_YOUSO[25]+1)+1];
  ****************************************************/

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(RF[i]);
  }
  free(RF);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(dRF[i]);
  }
  free(dRF);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(AF[i]);
  }
  free(AF);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(dAFQ[i]);
  }
  free(dAFQ);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(dAFP[i]);
  }
  free(dAFP);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(ddAFPP[i]);
  }
  free(ddAFPP);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(ddAFQQ[i]);
  }
  free(ddAFQQ);

  for (i=0; i<(List_YOUSO[25]+1); i++){
    free(ddAFQP[i]);
  }
  free(ddAFQP);

}



