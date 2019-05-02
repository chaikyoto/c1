//偏微分方程式の陽的解法
//中心差分法
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);
double writefile(double *u, double t, double Nx, double dx, int num, int inv);

int main(void){

  double xmax, tmax, dx, dt;
  int Nx, Nt,i,n,x,t,j,num,inv;
  double *u;

  //値域
  xmax = 1.0;
  tmax = 4.0e-2;

  //刻み幅
  dx = 0.01;
  dt = 1.0e-8;

  //格子
  Nx = xmax/dx+1;
  Nt = tmax/dt+1;

  //書き出し間隔
  inv = 20000;

  //あるtにおけるやつを計算するu(x)
  u = malloc(sizeof(double *) * Nx);

  //t=0の初期条件
  t = 0;
  for(i=0;i<Nx;i++){
    u[i] = f(i*dx);
  }

  //初期値を書き出し
  num = 0;
  writefile(u, t, Nx, dx, num, inv);

  //時間を進める
  for(j=0;j<Nt;j++){
    t += dt;
    for(i=1;i<Nx;i++){
      u[i] += ((i+1)*dx*u[i+1]-2*i*dx*u[i]+(i-1)*dx*u[i-1])*dt/(pow(i*dx,3)*dx*dx);
    }
    if(j%inv == 0){
      num++;
      writefile(u, t, Nx, dx, num, inv);
    }
  }

  return 0;
}

double f(double x){
  //return exp(-(x-0.5)*(x-0.5)/0.25);
  if(fabs(x-0.5)<=0.1){
    return 1.0;
  }else{
    return 0.0;
  }
}

//ファイル書き出し関数
double writefile(double *u, double t, double Nx, double dx, int num, int inv){
  //書き出し
  FILE *fp;
  char fname[50];
  int k;

  sprintf(fname, "dat/accd_n%05d.dat", num);
  fp = fopen(fname, "w");

  if(fp == NULL){
    printf("error\n");
    exit(1);
  }
  //fprintf(fp, "t=%12.4e inv=%03d\n", t, inv);


  for(k=0;k<Nx;k++){
    fprintf(fp, "%12.4e %12.4e\n", k*dx, u[k]);
  }
  fclose(fp);

  return 0;
}
