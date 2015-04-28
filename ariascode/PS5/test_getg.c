#define Nmx 40000
#define Itmx 100
#define alpha 0.25
#include "nrutil.h"
#include "p480.h"
#include <stdio.h>
#include <math.h>
#include "physics.h"

main()
{
  /* Change of variables info */
  double *r,*dr;
  int N;

  /* Physics variables */
  double *Rho,*g,Z;
  int lmax,*nmax,nmaxmax;
  double **F;
  double Etot;

  /* Working variables */
  int n,l,k;
  double x;

  /* Solver iteration variables */
  int it;
  /* Value of pi */
  const double pi=4.*atan(1.);

  /* Specs for O */
  Z=8.;
  lmax=1;

  nmax=ivector(0,lmax);
  nmax[0]=1;
  nmax[1]=0;

  nmaxmax=0;
  for (l=0; l<=lmax; l++)
    if (nmax[l]>nmaxmax) nmaxmax=nmax[l];

  F=dmatrix(0,lmax,0,nmaxmax);
  F[0][0]=2.;
  F[0][1]=2.;
  F[1][0]=4.;
  
  /* The rest is now general for ANY case */
  Rho=dvector(0,Nmx);
  g=dvector(0,Nmx);

  /* Grid vectors */
  r=dvector(0,Nmx); 
  dr=dvector(0,Nmx);

  N=4000;
  /* Set up grid */
  for (k=0; k<=N; k++) {
    x=((double) k)/((double) N);
    r[k]=1/(1-x)-1-x-x*x-x*x*x;
    dr[k]=1/(1-x)/(1-x)-1-2*x-3*x*x;
  }
  dr[N]=0.;

  /* Initialize charge density */
  for (k=0; k<=N; k++) 
    Rho[k]=0.;
  
  /* Iteration loop */
  for (it=1; it<=Itmx; it++) {
    Etot=getg(g,Rho,Z,lmax,nmax,nmaxmax,F,r,dr,N);
    printf("%4d  %20.12f   %f\n",
           it,fabs(Etot-(-74.473076803203738)), log(fabs(Etot-(-74.473076803203738))));

    for (k=0; k<=N; k++)
      Rho[k]=Rho[k]+alpha*g[k];
  }

  /* Be a good citizen and clean up... */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);

  free_dmatrix(F,0,lmax,0,nmaxmax);
  free_dvector(Rho,0,Nmx);
  free_dvector(g,0,Nmx);
}
