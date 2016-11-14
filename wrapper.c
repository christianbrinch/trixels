#include <stdio.h>
#include <stdlib.h>

int main(){
  extern void fourier_();

  int n,t,d, *tr, i,j;
  float *x, *y, *z;
  FILE *fp;

  if((fp=fopen("temp_out", "r")) == NULL) {
     printf("Cannot open file.\n");
	 exit(1);
  }

  fscanf(fp, "%d", &n);
  fscanf(fp, "%d", &t);
  fscanf(fp, "%d", &d);

  x=malloc(sizeof(float)*n);
  y=malloc(sizeof(float)*n);
  z=malloc(sizeof(float)*n);
  tr=malloc(sizeof(int)*3*t);

  for(i=0;i<n;i++){
	fscanf(fp, "%f", &x[i]);
  }
  for(i=0;i<n;i++){
	fscanf(fp, "%f", &y[i]);
  }
  for(i=0;i<n;i++){
	fscanf(fp, "%f", &z[i]);
  }
  for(i=0;i<t;i++){
	for(j=0;j<3;j++) fscanf(fp, "%d", &tr[j+3*i]);
  }

  

  fourier_(&n,&t,&d,x,y,z,tr);

  free(x);
  free(y);
  free(z);
  free(tr);

  return 1;
}
