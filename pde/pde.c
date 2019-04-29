#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);

int main(void){

  double xmax, tmax, dx, dt;
  int Nx, Nt, i,n,x,t,j;
  double *u;

  //値域
  xmax = 1.0;
  tmax = 1.0;
  
  //刻み幅
  dx = 0.01;
  dt = 0.01;

  //格子
  Nx = xmax/dx+1;
  Nt = tmax/dt+1;

  u = malloc(sizeof(double *) * Nx * Nt);

  for(n=0;n<Nt;n++){
    u[0+n] = f(0.0);
    u[Nt*(Nx-1)+n] = f(xmax);
  }

  for(x=1;x<(Nx-1);x++){
    for(t=0;t<(Nt-1);t++){
      u[x*Nt + t+1] = u[x*Nt + t] + t*(u[(x+1)*Nt + t]-2*u[x*Nt + t] + u[(x-1)*Nt + t])/(dx*dx);
    }
  }

  for(i=0;i<Nt;i++){
    for(j=0;j<Nx;j++){
      printf("%12.4e",u[i*Nt + j]);
    }
    printf("\n");
  }
  

  return 0;
}

double f(double x){
  return exp(-(x-0.5)*(x-0.5)/0.25);
}
