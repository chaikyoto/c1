#include<stdio.h>
#include<math.h>
#include<stdlib.h>

double f (int i, double x, double *y );

void runge4(double x, double *y, int n, double h, double funcpt(int i, double x, double *y));

int main(void)
{
  double x, y[8], h, xmax;
  int n, imax, i;

  x = 0.0;
  xmax = 1000.0;
  y[0] = 1.0;
  y[1] = 0.0;
  y[2] = -1.0;
  y[3] = 0.0;
  y[4] = 0.0;
  y[5] = 1.0;
  y[6] = 0.0;
  y[7] = -1.0;

  n = 8;
  h = 0.1;
  imax = xmax / h + 1;

  printf("%12.4e %12.4e %12.4e %12.4e %12.4e\n", x, y[0], y[1],y[2],y[3]);

  for(i=1; i<=imax; i++)
    {
      runge4(x, y, n, h, f);

      x = x + h;

      printf("%12.4e %12.4e %12.4e %12.4e %12.4e\n", x, y[0], y[1],y[2],y[3]);

    }
  return 0;
}

double f( int i, double x, double *y)
{
  double eps = 0.1;
  double R1 = pow(sqrt(y[0]*y[0] + y[1]*y[1] + eps*eps),3);
  double R2 = pow(sqrt(y[2]*y[2] + y[3]*y[3] + eps*eps),3);
  double R12 = pow(sqrt((y[2]-y[0])*(y[2]-y[0]) + (y[3]-y[1])*(y[3]-y[1]) + eps*eps),3);

  if (i == 0)
    return( y[4] );
  else if (i == 1 )
    return( y[5] );
  else if (i == 2)
    return( y[6] );
  else if (i == 3)
    return( y[7] );
  else if (i == 4)
    return( -(1/R1 + 1/R12)*y[0] + 1/R12 * y[2] );
  else if (i == 5)
    return( -(1/R1 + 1/R12)*y[1] + 1/R12 * y[3] );
  else if (i == 6)
    return( -(1/R2 + 1/R12)*y[1] + 1/R12 * y[3] );
  else if (i == 7)
    return( -(1/R2 + 1/R12)*y[0] + 1/R12 * y[2] );
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
