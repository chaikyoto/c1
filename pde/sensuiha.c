//偏微分方程式の陽的解法
//中心差分法
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double f(double x);
double writefile(double *h, double t, double Nx, double dx, int num, int inv);

int main(void){

  double xmax, tmax, dx, dt;
  int Nx, Nt,i,n,x,t,j,num,inv,k,l;
  double *u, *h, H, g;
  double *tmpu, *tmph, *b1, *b2;
  double **A, *a;


  //値域
  xmax = 1.0e4;
  tmax = 100;

  //刻み幅
  dx = 100;
  dt = 1.0;

  //格子
  Nx = xmax/dx+1;
  Nt = tmax/dt+1;

  //書き出し間隔
  inv = 1;

  //あるtにおけるやつを計算するので、u(x)
  u = malloc(sizeof(double *) * Nx);
  h = malloc(sizeof(double *) * Nx);
  tmpu = malloc(sizeof(double *) * Nx);
  tmph = malloc(sizeof(double *) * Nx);
  b1 = malloc(sizeof(double *) * Nx);
  b2 = malloc(sizeof(double *) * Nx);

  //t=0の初期条件
  t = 0;
  for(i=0;i<Nx;i++){
    h[i] = f(i*dx);
    u[i] = 0.0;
  }

  //深さH
  H = 1000;
  //重力加速度
  g = 9.8;

  //初期値を書き出し
  num = 0;
  writefile(h, t, Nx, dx, num, inv);

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

  //行列の中身に書き込み
  for(i=1;i<Nx-1;i++){
      A[i][i-1] = -1.0;
      A[i][i] = 2.0;
      A[i][i+1] = -1.0;
  }

  free(a);

  //時間を進める
  for(j=0;j<Nt;j++){

    for(k=0;k<Nt;k++){
      tmpu[k] = u[k];
      tmph[k] = h[k];
    }

    t += dt;

    u[0] = u[Nx-1];
    h[0] = h[Nx-1];

    for(i=0;i<Nx;i++){
      b1[i] = u[i] - tmpu[i]+0.1;
      b2[i] = h[i] - tmph[i]+0.1;
    }

    for(i=0;i<Nx;i++){
      SOR(A, u, b2, -dx*dx/H/dt, Nx);
    }
    for(i=0;i<Nx;i++){
      SOR(A, h, b1, -dx*dx/g/dt, Nx);
    }


    if(j%inv == 0){
      num++;
      writefile(h, t, Nx, dx, num, inv);
    }
  }

  free(tmpu);
  free(tmph);
  free(b1);
  free(b2);
  free(A);

  return 0;
}

double f(double x){
  //return exp(-(x-0.5)*(x-0.5)/0.25);

  return pow(sin(3.1415*x/10000),100);

}

//SOR法
void SOR(double **A, double *x, double *b, double fac, int N){

  int i, j, k, itr;
  double omega, eps, xtilde, norm, *newx;

  newx = malloc(sizeof(double *) * N);

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
      xtilde = b[i];
      for(j=0;j<i;j++){
        xtilde -= A[i][j]*newx[j];
      }
      for(k=i+1;k<N;k++){
        xtilde -= A[i][k]*newx[k];
      }
      xtilde = xtilde / A[i][i];
      newx[i] += omega * (xtilde - newx[i]);

      //収束判定用
      norm = norm + fabs(xtilde-newx[i]);
    }

    if(itr>1.0e5){
      printf("ERROR\n");
      exit(1);
    }

  }while(norm>eps);

  //値を更新
  for(i=0;i<N;i++){
    x[i] = fac*newx[i];
  }

  free(newx);

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
