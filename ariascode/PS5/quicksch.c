#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

void derivs_Schrodinger(double x, double y[], double dydx[], int k, double T[], double r[], double dr[], int N){
  dydx[1] = y[2]*dr[k];
  dydx[2] = -2.*T[k]*y[1]*dr[k];
}

#define Nmx 200

main()
{
  /* Grid vectors */
  double *r,*dr,*T;

  /* Working variables */
  int i,N;
  double u;
  double *y,*dydx; /* NR vector for diff eq's */
  double x,h;
  int k,dk;

  double E;

  /* Value of pi */
  double pi;
  pi=4.*atan(1.);

  /* Allocate NR vectors for maximum size used */
  r=dvector(0,Nmx);
  dr=dvector(0,Nmx);
  T=dvector(0,Nmx);
  y=dvector(1,2);
  dydx=dvector(1,2);

  /* Loop over different integration mesh sizes */
  N=Nmx;

  while (1) {
    printf("E: "); scanf("%lf",&E);

    /* Construct change of variables information and kinitic energy
       for solving for phi */
    for (i=0; i<=N; i++) {
      u=((double) i)/((double) N);
      r[i]=1./(1.-u)-1.;
      dr[i]=1./(1.-u)/(1.-u);
      T[i]=(E+1./r[i]);
    }

    /* Runga-Kutta solution using rk4p480(). */
    /* Set up step sizes (h,dk) and initial conditions here ... */
    dk = -2;
    h = -2./N;
    y[1]=0;
    y[2]=1e-20;
    //    r[N]=0.;
    dr[N]=0.;
    //T[0]=0.;

    for (k=N; k>=2; k+=dk) {
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
    
  }
    printf("N=%d; E, Psi[0]: %20.16f %e\n",N,E,y[1]);
 
  }
  /* Deallocate NR vectors: Should always clean up space at the end. */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dvector(T,0,Nmx);
  free_dvector(y,1,1);
  free_dvector(dydx,1,1);
  }


