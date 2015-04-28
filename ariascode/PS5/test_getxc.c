#include "p480.h"
#include "nrutil.h"
#define Nmx 40000
# include <math.h>

main()
{
   /* Change of variables info */
  double *r,*dr;
  int N;

  /* Physics variables */
  double *V,*Rho,*Rhonew,*phi,*Vxc,*Depsxc,*integrand,Z;
  int lmax,*nmax,nmaxmax;
  double **E,***Psi,**F;
  double Etot;

  /* Working variables */
  int n,l,k;
  double x;

  /* Solver iteration varibles */
  int it;

  /* Value of pi */
  const double pi=4.*atan(1.);

  /* Specs for hydrogen */
  Z=1.;
  lmax=0;

  nmax=ivector(0,lmax);
  nmax[0]=0;

  nmaxmax=0;
  for (l=0; l<=lmax; l++)
    if (nmax[l]>nmaxmax) nmaxmax=nmax[l];

  F=dmatrix(0,lmax,0,nmaxmax);
    F[0][0]=1.;

  /* The rest is now general for ANY case */
  E=dmatrix(0,lmax,0,nmaxmax); /* Make space for E's and Psi's */
  Psi=d3tensor(0,lmax,0,nmaxmax,0,Nmx);
  Rho=dvector(0,Nmx);
  Rhonew=dvector(0,Nmx);
  phi=dvector(0,Nmx);
  Vxc=dvector(0,Nmx);
  Depsxc=dvector(0,Nmx);

  integrand=dvector(0,Nmx);

  /* Grid vectors */
  r=dvector(0,Nmx); 
  dr=dvector(0,Nmx);
  V=dvector(0,Nmx);

  for (N=40; N<=Nmx; N*=10) {

    /* Set up grid */
    for (k=0; k<=N; k++) {
      x=((double) k)/((double) N);
      r[k]=1/(1-x)-1-x;
      dr[k]=1/(1-x)/(1-x)-1;
    }
    dr[N]=0.;

    /* Initialize charge density */
    for (k=0; k<=N; k++)
      Rho[k]=0.;
    
    /* Test section */
    {
      /* Make potential from zero charge (debugs NaNs, etc.) */
      getphi(phi,Rho,r,dr,N);
      getVxc(Vxc,Rho,r,dr,N);
      for (k=1; k<=N; k++)
        V[k]=-Z/r[k]+phi[k]+Vxc[k];
      V[0]=0.;
      
      /* Get H 1s wave function, and DENSITY */
      getallEs(E,lmax,nmax,Z,V,r,dr,N);
      getallPsis(Psi,E,lmax,nmax,V,r,dr,N);
      getRho(Rhonew,Psi,F,lmax,nmax,N);

      /* Compute Exc energy correction */
      getDepsxc(Depsxc,Rhonew,r,dr,N);
      getVxc(Vxc,Rhonew,r,dr,N);

      for (k=0; k<=N; k++)
        integrand[k]=Vxc[k]*Rhonew[k];
      printf("N=%6d <Vxc>=%15.10f ",N,
             simpint(integrand,r,dr,N)
             );
      for (k=0; k<=N; k++)
        integrand[k]=(Depsxc[k])*Rhonew[k];
      printf("<DExc>=%15.10f\n",
             simpint(integrand,r,dr,N)
             );
    }
  }

  /* Be a good citizen and clean up... */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dvector(V,0,Nmx);

  free_dmatrix(E,0,lmax,0,nmaxmax);
  free_d3tensor(Psi,0,lmax,0,nmaxmax,0,Nmx);
  free_dvector(Rho,0,Nmx);
  free_dvector(Rhonew,0,Nmx);
  free_dvector(phi,0,Nmx);
  free_dvector(Vxc,0,Nmx);
  free_dvector(Depsxc,0,Nmx);
}

