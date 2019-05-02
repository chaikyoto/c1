//SOR法で連立一次方程式Ax=bを解く
//固有値がすべて正の対称行列じゃないと収束しない！
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main(void){

  double **A, *a, *b, *x;
  int i, j, k, t, N;
  double omega, eps, xtilde, norm;

  N = 3;

  //メモリの動的確保
  A = malloc(sizeof(double *) * N);
  a = malloc(sizeof(double) * N * N);
  b = malloc(sizeof(double *) * N);
  x = malloc(sizeof(double *) * N);

  for(i=0;i<N;i++){
    A[i] = a + i * N;
  }

  A[0][0] = 1.0;
  A[0][1] = 2.0;
  A[0][2] = 3.0;

  A[1][0] = 1.0;
  A[1][1] = 3.0;
  A[1][2] = 4.0;

  A[2][0] = 2.0;
  A[2][1] = 3.0;
  A[2][2] = 8.0;

  b[0] = 3.0;
  b[1] = 2.0;
  b[2] = 1.0;


  //行列の表示
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      printf("%e  ", A[i][j]);
      if(j==N-1){
        printf("\n");
      }
    }
  }

  //反復回数t
  t = 0;

  //加速パラメータ
  omega = 1.5;

  //収束判定
  eps = 1.0e-8;

  do{
    t++;
    norm = 0;

    for(i=0;i<N;i++){
      xtilde = b[i];
      for(j=0;j<i;j++){
        xtilde -= A[i][j]*x[j];
      }
      for(k=i+1;k<N;k++){
        xtilde -= A[i][k]*x[k];
      }
      xtilde = xtilde / A[i][i];
      x[i] += omega * (xtilde - x[i]);

      //収束判定用
      norm = norm + fabs(xtilde-x[i]);
    }

  }while(norm>eps);

  //解の表示
  for(i=0;i<N;++i){
    printf("%12.4e\n", x[i]);
  }


  free(a);
  free(A);

  return 0;
}
