
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



int main(int argc, char *argv[]) 
{ 
  char filename[100];
  char ctmp[100],ctmp1[100];
  char dmyss[256];
  int Ngrid1,Ngrid2,Ngrid3;
  FILE *fp;
  double x,y,z;
  double val1,val2,val3;
  int i,j,k;


  strcpy(filename,argv[1]);
  //printf("%s \n",filename);

  if ((fp = fopen(filename,"r")) == NULL){
    printf("cannot open %s \n", filename);
    exit(0);
  }

  /* first line */
  fscanf(fp,"%s %s %d %d %d",ctmp,ctmp1,&Ngrid1,&Ngrid2,&Ngrid3);
  fgets(dmyss,256,fp);
  //printf("%d %d %d \n",Ngrid1,Ngrid2,Ngrid3);
  //
  /* second line */
  fscanf(fp,"%s",ctmp);
  fgets(dmyss,256,fp);

  for (i=0;i<Ngrid1;i++){
    for (j=0;j<Ngrid2;j++){
      for (k=0;k<Ngrid3;k++){
        fscanf(fp," %lf %lf %lf %lf",&x,&y,&z,&val1);
        //fscanf(fp,"%lf %lf %lf %lf %lf %lf",&x,&y,&z,&val1,&val2,&val3);
        fgets(dmyss,256,fp);
        //printf("%lf %lf %lf %lf \n",x,y,z,val1);
        //printf("%14.4E %14.4E %14.4E %14.4E \n",x,y,z,val1);
        printf("%14.4E",val1);
        if ((k+1)%6==0) { printf("\n"); }
      }
      /* avoid double \n\n when Ngrid3%6 == 0  */
      if (Ngrid3%6!=0) { printf("\n"); }
    }
  }




  fclose(fp);

}
