 //ものすごく単純なN体シミュレーション
//imax個のファイルを作ります
//vimでうってみたよ。テストテスト。
//test2
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define N 100 //粒子数

void runge4(double t, double *y, double *m, int n, double h, double funcpt(int i, double t, double *y, double *m));
double f(int i, double t, double *y, double *m);
double g(int k, double *x, double *m, int p);
//0から1までの乱数を作る関数
double ran();
double plamer();

int main(void){
  srand(time(NULL));

  //ファイル構造体
  FILE *fp;

  //x:position, v:velocity, m:mass
  double t, m[N],y[4*N], h, tmax, xmax, ymax, rmax, R, theta;
  int n, imax, i, j, l, k, q, interval;
  char fname[50];

  t = 0.0;
  tmax = 1.0;

  h = 0.01;

  imax = tmax / h;

  //xmax = 10.0;
  //ymax = 10.0;

  rmax = 1.0;

  //書き出し間隔
  interval = 1;

/*
  m[0] = 10000.0;
  y[0] = 0.0;
  y[N] = 0.0;
  y[2*N] = 0.0;
  y[3*N] = 0.0;
  */

  //initial conditions
  for(i=0; i<N; i++){
    //R = plamer();
    m[i] = 0.1;
    R = ran()*rmax;
    theta = ran()*2*M_PI;
    y[i] = R*cos(theta);
    y[N+i] = R*sin(theta);
    y[2*N+i] = -sin(theta)/sqrt(R);
    y[3*N+i] = cos(theta)/sqrt(R);
    //y[2*N+i] = 0.0;
    //y[3*N+i] = 0.0;
  }


  //ファイル書き込み
  fp = fopen("nbody_00000.dat","w");
  if(fp == NULL){
    printf("ERROR\n");
    exit(1);
  }
  for(i=0; i<N; i++){
    fprintf(fp,"%12.4e %12.4e\n",y[i],y[N+i]);
  }
  fclose(fp);

  //printf("%12.4e %12.4e %12.4e %12.4e\n",y[0],y[2],y[1],y[3]);

  //計算ループ
  int filenum = 0;

  for(i=1; i<=imax; i++){

    runge4(t, y, m, 4*N, h, f);

    t = t + h;

    //ファイル書き込み
    if(i%interval==0){
      filenum++;
      sprintf(fname,"nbody_%05d.dat",filenum);
      fp = fopen(fname,"w");
      if(fp == NULL){
        printf("ERROR\n");
        exit(1);
      }
      for(j=0; j<N; j++){
        fprintf(fp, "%12.4e %12.4e\n",y[j],y[N+j]);
      }
      fclose(fp);
    }
    //printf("%12.4e %12.4e %12.4e %12.4e\n",y[0],y[2],y[1],y[3]);

  }

}

double plamer(){
  double r0 = 1;
  double xp,yp;
  while(1){
    xp = 3*ran();
    yp = 0.25*ran();
    if(yp<0.25*r0*r0/pow((xp*xp+r0*r0),2.5)){
      break;
    }
  }
  return xp;
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
  double eps = 1.0e-3;
  double sum, r2;
  int j;

  sum = 0;

  for(j=0; j<N; j++){
    if(j != k){
      r2 = pow((y[j]-y[k]), 2.0) + pow((y[N+j]-y[N+k]), 2.0);
      //printf("r2: %12.4e\n",r2);
      sum += m[j]*(y[N*p+j] - y[N*p+k])/pow(r2+eps*eps, 1.5);
    }
  }
  //printf("sum: %12.4e\n",sum);
  return sum;
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
