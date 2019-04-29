//二体問題のシミュレーション
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double f (int i, double x, double *y, double eps);

void runge4(double x, double *y, int n, double h, double epsilon, double funcpt(int i, double x, double *y, double eps));

int main(void)
{
  double x, y[4], h, xmax, epsilon;
  int n, imax, i;

  x = 0.0;
  xmax = 10.0;
  y[0] = 1.0;
  y[1] = 0.0;
  y[2] = -1.0;
  y[3] = 2.0;

  n = 4;
  h = 0.1;
  epsilon = 1.0e-3;
  //epsilon = 0.0;
  imax = xmax / h + 1;

  printf("%12.4e %12.4e %12.4e\n", x, y[0], y[1]);

  for(i=1; i<=imax; i++)
    {
      runge4(x, y, n, h, epsilon, f);

      x = x + h;

  printf("%12.4e %12.4e %12.4e\n", x, y[0], y[1]);

    }
  return 0;
}

double f( int i, double x, double *y, double eps)
{
  double r = sqrt(y[0]*y[0] + y[1]*y[1] + eps*eps);
  if (i == 0)
    return( y[2] );
  else if (i == 1 )
    return( y[3] );
  else if (i == 2)
    return( -y[0]/(r*r*r) );
  else if (i == 3)
    return( -y[1]/(r*r*r) );
  return 0;
}

void runge4(double x, double *y, int n, double h, double epsilon, double funcpt(int i, double x, double *y, double eps))
{
  double eps = epsilon;
  double *k1pt, *k2pt, *k3pt, *k4pt, *work;
  int i;

  k1pt = (double *)malloc(n*8);
  k2pt = (double *)malloc(n*8);
  k3pt = (double *)malloc(n*8);
  k4pt = (double *)malloc(n*8);
  work = (double *)malloc(n*8);

  for(i = 0; i < n; i++)
    k1pt[i] = h* funcpt(i,x,y, eps);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k1pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k2pt[i] = h * funcpt( i, x+h/2.0, work, eps);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k2pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k3pt[i] = h * funcpt( i, x+h/2.0, work, eps);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k3pt[i];
  for(i = 0; i < n; i++)
    k4pt[i] = h * funcpt( i, x+h, work, eps);
  for(i = 0; i < n; i++ )
    y[i] += ( k1pt[i] + 2.0 * k2pt[i] + 2.0 * k3pt[i] + k4pt[i] ) / 6.0;

  free(k1pt);
  free(k2pt);
  free(k3pt);
  free(k4pt);
  free(work);
}
