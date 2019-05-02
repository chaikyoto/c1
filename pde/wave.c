//浅水波方程式
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double writefile(double *h, double t, double Nx, double dx, int num, int inv);
void rk4(double *x, double t, double dt, double *f);

int main(void){

  for(i=0;i<Nt;i++){
    t += dt;
    rk4(i, u, t, dt, du);
    rk4(i, h, t, dt, dh);

    if(i%inv == 0){
      num++;
      writefile(h, t, Nx, dx, num, inv);
    }

  }

}

double du(double x, double t){
  return -g*(h[i+2]-h[i])/(2*dx);
}

void rk4(int i, double *x, double t, double h, double funcpt(double x, double t)){
  double k1, k2, k3, k4;
  k1 = h*f(x[i],t);
  k2 = h*f(x[i] + 0.5*k1, t + 0.5*h);
  k3 = h*f(x[i] + 0.5*k2, t + 0.5*h);
  k4 = h*f(x[i] + k3, t + h);
  x[i] += (k1 + 2*k2 + 2*k3 + k4) / 6;
}

//ファイル書き出し関数
double writefile(double *h, double t, double Nx, double dx, int num, int inv){
  //書き出し
  FILE *fp;
  char fname[50];
  int k;

  sprintf(fname, "sensui/2_n%05d.dat", num);
  fp = fopen(fname, "w");

  if(fp == NULL){
    printf("error\n");
    exit(1);
  }
  //fprintf(fp, "t=%12.4e inv=%03d\n", t, inv);


  for(k=0;k<Nx;k++){
    fprintf(fp, "%12.4e %12.4e\n", k*dx, h[k]);
  }
  fclose(fp);

  return 0;
}
