#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

void derivs_poisson(double x, double y[], double dydx[], int k, double rho[], double r[], double dr[], int N){
  
  double pi;
  pi = 4.*atan(1.);
  dydx[1] = y[2]*dr[k];
  dydx[2] = -4.*pi*rho[k]*r[k]*dr[k];
  //  printf("dydx1dydx2: %20.16f %20.16f\n",dydx[1], dydx[2]);
}
#define Nmx 200000
main()
{
  /* Grid vectors */
  double *r,*dr,*rho;
  
  /* Working variables */
  int i,N;
  double u;
  double *y,*dydx; /* NR vector for diff eq's */
  double x,h;
  int k,dk;
  
  double rms;
  
  /* Value of pi */
  double pi;
  pi=4.*atan(1.);
  
  /* Allocate NR vectors for maximum size used */
  r=dvector(0,Nmx);
  dr=dvector(0,Nmx);
  rho=dvector(0,Nmx);
  y=dvector(1,2);
  dydx=dvector(1,2);
  
  /* Loop over different integration mesh sizes */
  for (N=2; N<=Nmx; N*=10) {
  //N=20;
    /* Construct change of variables information and
       electron density "rho" for solving for phi */
    for (i=0; i<=N; i++) {
      u=((double) i)/((double) N);
      r[i]=1./(1.-u)-1.;
      dr[i]=1./(1.-u)/(1.-u);
      rho[i]=exp(-2.*r[i])/pi;
      //  printf("r,dr,rho : %20.16f %20.16f %20.16f\n",r[i],dr[i],rho[i]);
    }
    
    /* Runge-Kutta solution using rk4p480(). */
    /* Set up initial step sizes (h,dk) and initial conditions ... */
    h = -2./(double)N;
    dk = -2;
    y[1]=1.;
    y[2]=0.;
    rho[N] = 0.;
    dr[N] = 0.;
    r[N] = 0.;
    rms=0.;
    //    printf("N, phi[0], rms error: %6d %20.16f %20.16f\n",N,y[1],y[2]);
    for (k=N; k>=2; k+=dk) {
      x=-k*h;
      derivs_poisson(x,y,dydx,k,rho,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_poisson,k,dk,rho,r,dr,N);
      /* This compares your solution with the analytic solution.  This
	 comparison must be made *after* you make the step with rk4p480.
	 
	 Replace the ...'s with the appropriate analytic expression
	 from Problem Set 1 in terms of r[k]!  Remember that, here, we
	 solve for Phi(r)=r*phi(r)
      */
      rms = rms + pow((y[1]-(1.-exp(-2.*r[k-2])-r[k-2]*exp(-2.*r[k-2])) ), 2.);
      // printf("N, phi[0], rms error: %6d %20.16f %20.16f\n",N,y[1],rms);
    }
    rms=sqrt(rms/(N/2));
    printf("N, phi[0], rms error: %6d %20.16f %20.16f\n",N,y[1],rms);
  
   }  
  /* Deallocate NR vectors: Should always clean up space at the end. */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dvector(rho,0,Nmx);
  free_dvector(y,1,1);
  free_dvector(dydx,1,1);
}

  
