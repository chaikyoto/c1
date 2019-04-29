#include<stdio.h>
#include<stdlib.h>

double f (int i, double x, double *y );

void runge4(double x, double *y, int n, double h, double funcpt(int i, double x, double *y));

int main(void)
{
  double x, y[2], h;
  int n, imax, i;

  x = 0.0;
  y[0] = 0.0;
  y[1] = 1.0;

  n = 2;
  h = 0.1;
  imax = 2.0 / h + 1;

  printf("%12.4e %12.4e%12.4e %12.4e%12.4e\n", x, y[0], x-0.5*x*x, y[1], 1.0-x);

  for(i=1; i<=imax; i++)
    {
      runge4(x, y, n, h, f);

      x = x + h;

      printf("%12.4e %12.4e%12.4e %12.4e%12.4e\n", x, y[0], x-0.5*x*x,y[1], 1.0-x);
    }
  return 0;
}

double f( int i, double x, double *y)
{
  if (i == 0 )
    return( y[1] );
  else if (i == 1 )
    return( -1.0 );
}

void runge4(double x, double *y, int n, double h, double funcpt(int i, double x, double *y))
{

  double *k1pt, *k2pt, *k3pt, *k4pt, *work;
  int i;

  k1pt = (double *)malloc(n*8);
  k2pt = (double *)malloc(n*8);
  k3pt = (double *)malloc(n*8);
  k4pt = (double *)malloc(n*8);
  work = (double *)malloc(n*8);

  for(i = 0; i < n; i++)
    k1pt[i] = h* funcpt(i,x,y);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k1pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k2pt[i] = h * funcpt( i, x+h/2.0, work );
  for(i = 0; i < n; i++)
    work[i] = y[i] + k2pt[i] / 2.0;
  for(i = 0; i < n; i++)
    k3pt[i] = h * funcpt( i, x+h/2.0, work);
  for(i = 0; i < n; i++)
    work[i] = y[i] + k3pt[i];
  for(i = 0; i < n; i++)
    k4pt[i] = h * funcpt( i, x+h, work );
  for(i = 0; i < n; i++ )
    y[i] += ( k1pt[i] + 2.0 * k2pt[i] + 2.0 * k3pt[i] + k4pt[i] ) / 6.0;

  free(k1pt);
  free(k2pt);
  free(k3pt);
  free(k4pt);
  free(work);
}
