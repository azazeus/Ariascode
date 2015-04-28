#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

#define Nmx 400
main()
{
  /* Change of variables info */
  double *r,*dr;
  int N;

  /* Physics variables */
  double *V,*Rho,Z;
  int lmax,*nmax,nmaxmax;
  double **E,***Psi,**F;

  /* Working variables */
  int n,l,k;
  double x;

  /* Value of pi */
  const double pi=4.*atan(1.);

  /* Specifications for carbon */
  Z=6.;
  lmax=1;

  nmax=ivector(0,lmax);
  nmax[0]=1;
  nmax[1]=0;

  nmaxmax=0; /* Find max of all nmax's */
  for (l=0; l<=lmax; l++)
    if (nmax[l]>nmaxmax) nmaxmax=nmax[l];

  F=dmatrix(0,lmax,0,nmaxmax);
  F[0][0]=2.; /* 2 electrons in 1s */
  F[0][1]=2.; /* 2 electrons in 2s */
  F[1][0]=2.; /* 2 electrons in 2p */

  /* The rest is now general for ANY case */
  E=dmatrix(0,lmax,0,nmaxmax); /* Make space for E's and Psi's */
  Psi=d3tensor(0,lmax,0,nmaxmax,0,Nmx);
  Rho=dvector(0,Nmx);

  /* Grid vectors */
  r=dvector(0,Nmx); 
  dr=dvector(0,Nmx);
  V=dvector(0,Nmx);

  /* Set up grid */
  N=400;
  for (k=0; k<=N; k++) {
    x=((double) k)/((double) N);
    r[k]=1./(1.-x)-1.-x-x*x;
    dr[k]=1./(1.-x)/(1.-x)-1.-2*x;
    V[k]=-Z/r[k];
  }
  V[0]=0.; dr[N]=0.;

  /* Test section */  
  getallEs(E,lmax,nmax,Z,V,r,dr,N);
  getallPsis(Psi,E,lmax,nmax,V,r,dr,N);
  getRho(Rho,Psi,F,lmax,nmax,N);
  printf("Total charge is: %20.15f\n",simpint(Rho,r,dr,N));

  /* Be a good citizen and clean up... */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dvector(V,0,Nmx);

  free_dmatrix(E,0,lmax,0,nmaxmax);
  free_dmatrix(F,0,lmax,0,nmaxmax);
  free_d3tensor(Psi,0,lmax,0,nmaxmax,0,Nmx);
  free_dvector(Rho,0,Nmx);
}
