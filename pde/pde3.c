//偏微分方程式の陰的解法
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);
void SOR(double **A, double *u, double beta, int N);
void writefile(double *u, double t, double Nx, double dx, int num, int inv);

int main(void){

  double xmax, tmax, dx, dt;
  int Nx, Nt,i,n,x,t,j,num,inv;
  double **A, *a, *b, *re;
  double beta, alpha;
  double *u;

  //値域
  xmax = 1.0;
  tmax = 1.0;

  //刻み幅
  dx = 0.01;
  dt = 1.0e-3;

  //格子
  Nx = xmax/dx+1;
  Nt = tmax/dt+1;

  //あるtにおけるやつを計算する。u(x)
  u = malloc(sizeof(double *) * Nx);

  //係数行列etcの確保
  A = malloc(sizeof(double *) * Nx);
  a = malloc(sizeof(double) * Nx * Nx);
  //二次にする
  for(i=0;i<Nx;i++){
    A[i] = a + i * Nx;
  }

  for(i=0;i<Nx;++i){
    for(j=0;j<Nx;++j){
      A[i][j]=0.0;
    }
  }

  //係数行列の中身
  beta = (dx*dx)/dt;
  alpha = beta + 2;

  //行列の中身に書き込み
  A[0][0] = beta;
  A[Nx-1][Nx-1] = beta;
  for(i=1;i<Nx-1;i++){
      A[i][i-1] = -1.0;
      A[i][i] = alpha;
      A[i][i+1] = -1.0;
  }

  //t=0の初期条件
  t = 0;
  for(i=0;i<Nx;i++){
    u[i] = f(i*dx);
  }

  //初期値を書き出し
  num = 0;
  writefile(u, t, Nx, dx, num, inv);

  //書き出し間隔
  inv = 1;

  //時間を進める
  for(j=0;j<Nt;j++){
    t += dt;

    SOR(A,u,beta,Nx);

    if(j%inv == 0){
      num++;
      writefile(u, t, Nx, dx, num, inv);
    }
  }

  free(a);
  free(A);
  free(u);

  return 0;
}

double f(double x){
  //return exp(-(x-0.5)*(x-0.5)/0.25);

  if(x<=0.5){
    return 1.0;
  }else{
    return 0.0;
  }

}

//SOR法
void SOR(double **A, double *u, double beta, int N){

  int i, j, k, itr;
  double omega, eps, utilde, norm, *newu;

  newu = malloc(sizeof(double *) * N);

  //反復回数itr
  itr = 0;

  //加速パラメータ
  omega = 1.5;

  //収束判定
  eps = 1.0e-10;

  do{
    itr++;
    norm = 0;

    for(i=0;i<N;i++){
      utilde = u[i];
      for(j=0;j<i;j++){
        utilde -= A[i][j]*newu[j];
      }
      for(k=i+1;k<N;k++){
        utilde -= A[i][k]*newu[k];
      }
      utilde = utilde / A[i][i];
      newu[i] += omega * (utilde - newu[i]);

      //収束判定用
      norm = norm + fabs(utilde-newu[i]);
    }

    if(itr>1.0e5){
      printf("ERROR\n");
      exit(1);
    }

  }while(norm>eps);

  //値を更新
  for(i=0;i<N;i++){
    u[i] = beta*newu[i];
  }

  free(newu);

}

//ファイル書き出し関数
void writefile(double *u, double t, double Nx, double dx, int num, int inv){
  //書き出し
  FILE *fp;
  char fname[50];
  int k;

  sprintf(fname, "dat/pde3test_n%05d.dat", num);
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

}
