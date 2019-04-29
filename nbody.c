//ものすごく単純なN体シミュレーション
//imax個のファイルを作ります
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <time.h>

#define N 3 //粒子数

void runge4(double t, double *y, double *m, int n, double h, double funcpt(int i, double t, double *y, double *m));
double f(int i, double t, double *y, double *m);
double g(int k, double *x, double *m, int p);
//0から1までの乱数を作る関数
double ran();

int main(void){
  srand(time(NULL));

  //ファイル構造体
  FILE *fp;

  //x:position, v:velocity, m:mass
  double t, x[N][2],v[N][2],m[N],y[4*N], h, tmax, xmax, ymax;
  int n, imax, i, j, l, k;
  char fname[50];

  t = 0.0;
  tmax = 100.0;

  h = 0.1;

  imax = tmax / h + 1;

  xmax = 10.0;
  ymax = 10.0;

  //initial conditions
  for(i=0; i<N; i++){
    m[i] = 1.0;
    v[i][0] = 1.0;
    v[i][1] = 0.0;
    x[i][0] = ran()*xmax;
    x[i][1] = ran()*ymax;
  }

/*
  fp = fopen("nbody_00000.dat","w");
  if(fp == NULL){
    printf("ERROR\n");
    exit(1);
  }
  for(i=0; i<N; i++){
    fprintf(fp,"%12.4e %12.4e\n",x[i][0],x[i][1]);
  }
  fclose(fp);
  */
  printf("%12.4e %12.4e\n",x[i][0],x[i][1]);

  //計算するために全部を１つの配列yに結合
  i=0;
  for(l=0; l<=1; l++){
    for(k=0; k<N; k++){
      y[i] = x[k][l];
      y[2*N+i] = v[k][l];
      i++;
    }
  }

  //計算ループ
  for(i=1; i<=imax; i++){

    runge4(t, y, m, 4*N, h, f);

    t = t + h;
/*
    //ファイル書き込み
    sprintf(fname,"nbody_%05d.dat",i);
    fp = fopen(fname,"w");
    if(fp == NULL){
      printf("ERROR\n");
      exit(1);
    }
    for(i=0; i<N; i++){
      fprintf(fp, "%12.4e %12.4e\n",y[i],y[N+i]);
    }
    fclose(fp);
*/
  printf("%12.4e %12.4e\n",y[i],y[N+i]);
  }

}

//0-1の乱数
double ran(){
  return (float)rand() / RAND_MAX;
}

//方程式
double f(int i, double t, double *y, double *m){

  //iが0からN-1のときv[0][0]からv[N-1][0]を、Nから2N-1のときv[0][1]からv[N-1][1]を
  if( i <= N-1 ){
    return y[2*N+i];
  }else if( i <= 2*N-1){
    return y[2*N+i];
  }else if( i <= 3*N-1){
    return g(i-2*N, y, m, 0);
  }else if( i <= 4*N-1){
    return g(i-3*N, y, m, 1);
  }

}

//右辺のサムネーション
double g(int k, double *y, double *m, int p){
  double eps = 0.1;
  double sum, r2;
  int j;

  for(j=0; j<N; j++){
    if(j != k){
      r2 = pow((y[j]-y[k]), 2.0) + pow((y[N+j]-y[N+k]), 2.0);
      sum += m[j]*(y[N*p+j] - y[N*p+k])/pow(r2+eps*eps, 1.5);
    }
  }
  return -1*sum;
}

void runge4(double t, double *y, double *m, int n, double h, double funcpt(int i, double t, double *y, double *m))
{

  double *k1pt, *k2pt, *k3pt, *k4pt, *work;
  int i;

  k1pt = (double *)malloc(n*8);
  k2pt = (double *)malloc(n*8);
  k3pt = (double *)malloc(n*8);
  k4pt = (double *)malloc(n*8);
  work = (double *)malloc(n*8);

  for(i = 0; i < n; i++)
    k1pt[i] = h* funcpt(i,t,y,m);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k1pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k2pt[i] = h * funcpt( i, t+h/2.0, work, m);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k2pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k3pt[i] = h * funcpt( i, t+h/2.0, work, m);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k3pt[i];
  for(i = 0; i < n; i++)
    k4pt[i] = h * funcpt( i, t+h, work, m);
  for(i = 0; i < n; i++ )
    y[i] += ( k1pt[i] + 2.0 * k2pt[i] + 2.0 * k3pt[i] + k4pt[i] ) / 6.0;

  free(k1pt);
  free(k2pt);
  free(k3pt);
  free(k4pt);
  free(work);
}
