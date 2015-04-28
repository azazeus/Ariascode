#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

#define Nmx 40000
main()
{
  /* Change of variables info */
  double *r,*dr;
  int N;

  /* Physics variables */
  double *V;
  double Z;
  int lmax;
  int nmaxmax,*nmax;
  double **E,*psi;

  int n,l;

  /* Working variables */ 
  double x,*integrand;
  int k;

  /* Value of pi */
  const double pi=4.*atan(1.);

  /* Specs for hydrogen */
  Z=1.; /* Charge on nucleus */
  lmax=0; /* Only up to 's' states (l=0) are filled in H */

  nmax=ivector(0,lmax); /* Declare space for max nodes for each l */
  nmax[0]=0; /* Specify that only up to the first (zero nodes) s state
                is filled in H */

  /* General loop for finding nmaxmax (needed to define arrays)...
     Note: This is a little silly for H, but is handy when you have
     some random element. */
  nmaxmax=0;
  for (l=0; l<=lmax; l++)
    if (nmax[l]>nmaxmax) nmaxmax=nmax[l];

  /* Allocate space for table for E's */
  E=dmatrix(0,lmax,0,nmaxmax);

  /* Convenient loop over N */
  for (N=40; N<=Nmx; N*=10) {
    /* Allocate space for objects with size which depends on N */
    psi=dvector(0,N); /* psi */
    r=dvector(0,N); /* Grid vectors */
    dr=dvector(0,N);
    V=dvector(0,N);
    
    integrand=dvector(0,N); /* Working grid for integrating things */

    /* Set up grid and -Z/r potential */
    for (k=0; k<=N; k++) {
      x=((double) k)/((double) N);
      r[k]=1./(1.-x)-1.-x-x*x-x*x*x;
      dr[k]=1./(1.-x)/(1.-x)-1.-2*x-3*x*x;
      V[k]=-Z/r[k];
    }

    /* Below gives proper limits at end points */
    V[0]=0.; dr[N]=0.;
    getallEs(E,lmax,nmax,Z,V,r,dr,N); /* Get E */

    //    printf("done with getallEs....");
    //printf("calling getPsi with l = %f\n",l);
    getPsi(psi,E[0][0],0,V,r,dr,N); /* Get Psi for E[l=0][n=0] */

    //    printf("done with getPsi....");

    for (k=0; k<=N-1; k++) /* Get integrand for root mean square (rms) error */
      integrand[k]=(psi[k]-2*r[k]*exp(-r[k]))*(psi[k]-2*r[k]*exp(-r[k]));

    printf("%6d %20.16f\n",N,sqrt(simpint(integrand,r,dr,N)));

    /* Free up space for N dependent objects before going to next N*/
    free_dvector(r,0,N);
    free_dvector(dr,0,N);
    free_dvector(V,0,N);
    
    free_dvector(psi,0,N);
    free_dvector(integrand,0,N);
  }

  /* Free up all remaining objects */
  free_dmatrix(E,0,lmax,0,nmaxmax);
  free_ivector(nmax,0,lmax); 
}
