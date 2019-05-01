//偏微分方程式の陽的解法
//中心差分法
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);

int main(void){

  double xmax, tmax, dx, dt;
  int Nx, Nt, i,n,x,t,j,I;
  double *u;

  //値域
  xmax = 1.0;
  tmax = 1.0;
  
  //刻み幅
  dx = 0.01;
  dt = 1.0e-5;

  //格子
  Nx = xmax/dx+1;
  Nt = tmax/dt+1;

  u = malloc(sizeof(double *) * Nx * Nt);

  

  //初期条件
  for(n=0;n<Nx;n++){
    u[Nt*n] = f(n);
  }

  
  //境界条件
  for(n=0;n<Nt;n++){
    u[0+n] = f(0.0);
    u[Nt*(Nx-1)+n] = f(xmax);
  }

  for(x=1;x<(Nx-1);x++){
    for(t=1;t<(Nt-1);t++){
      u[x*Nt + t+1] = u[x*Nt + t] + dt*(u[(x+1)*Nt + t]-2*u[x*Nt + t] + u[(x-1)*Nt + t])/(dx*dx);
    }
  }

  /**
  for(i=0;i<Nx;i++){
    for(j=0;j<Nt;j++){
      printf("%12.2e",u[i*Nt + j]);
    }
    printf("\n");
  }
  **/

  //書き出し
  FILE *fp;
  int filenum = 0;
  int inv = 10;
  char fname[50];

  for(j=0;j<Nt;j++){

    if(j%inv==0){
      
      sprintf(fname, "dat/pdetest_inv%03d_n%05d.dat", inv, filenum);
      fp = fopen(fname, "w");

      if(fp == NULL){
        printf("error\n");
        exit(1);
       }

      for(i=0;i<Nx;i++){
        fprintf(fp, "%12.4e %12.4e\n", i*dx, u[i*Nt + j]);
      }
      fclose(fp);
      filenum++;
    }
  }
  

  return 0;
}

double f(double x){
  return sin(x*2*3.141592);
  //return 1.0;
}
