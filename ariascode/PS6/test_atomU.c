#include "p480.h"
#include "nrutil.h"
#define Nmx 40000
#include <math.h>
#include <stdio.h>
#include "physics.h"

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
  const double alpha = 0.35;

  /* Working variables */
   int n,l,k;
  double x;
  
  /* Solver iteration varibles */
   int it=0;
  
  /* Value of pi */
   const double pi=4.*atan(1.);
  
  /* Specs for Sb */
  Z=92.;
  lmax=3;
         
  nmax=ivector(0,lmax);
  nmax[0]=6;
  nmax[1]=4;
  nmax[2]=3;
  nmax[3]=1;

  nmaxmax=0;
  for (l=0; l<=lmax; l++)
   if (nmax[l]>nmaxmax) nmaxmax=nmax[l];

    F=dmatrix(0,lmax,0,nmaxmax);
  F[0][0]=2.; /* 1s */
  F[0][1]=2.; /* 2s */
  F[0][2]=2.;
  F[0][3]=2.; /* 1s */
  F[0][4]=2.; /* 2s */
  F[0][5]=2.;
  F[0][6]=2.;


  F[1][0]=6.; /* 2p */
  F[1][1]=6.; /* 3p */
  F[1][2]=6.; /* 4p */
  F[1][3]=6.; /* 5p */
  F[1][4]=6.; /* 5p */

  F[2][0]=10.; /* 2p */
  F[2][1]=10.; /* 3p */
  F[2][2]=10.; /* 4p */
  F[2][3]=1.; /* 5p */

  F[3][0]=14.; /* 2p */
  F[3][1]=3.; /* 3p */


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
   while(1){
     it++;
     /* Make potential from zero charge (debugs NaNs, etc.) */
      getphi(phi,Rho,r,dr,N);
      getVxc(Vxc,Rho,r,dr,N);
      for (k=0; k<=N; k++)
	V[k]=-Z/r[k]+phi[k]+Vxc[k];
      V[0]=0.;
                                                                                
      /* Get H 1s wave function, and DENSITY */
      getallEs(E,lmax,nmax,Z,V,r,dr,N);
      getallPsis(Psi,E,lmax,nmax,V,r,dr,N);
      getRho(Rhonew,Psi,F,lmax,nmax,N);
     for(k=0;k<=N;k++)Rho[k]=(1. - alpha)*Rho[k]+alpha*Rhonew[k];

      /* Compute and output total energy */
      
      /* Get correction to sum of electron energies */
        getDepsxc(Depsxc,Rho,r,dr,N);
	for (k=0; k<=N; k++)
	  integrand[k]=(-0.5*phi[k]+Depsxc[k])*Rho[k];
	Etot=simpint(integrand,r,dr,N);

      /* Add on the sum of the electron energies times the occupancies */
	for (l=0; l<=lmax; l++)
	  for (n=0; n<=nmax[l]; n++)
	    Etot+=F[l][n]*E[l][n];
      
	printf("%d \t Etot: %20.15f\t E (from NIST): %f \t error: %f\n",it, Etot, -25658.417889, fabs(Etot+25658.417889));

   }
    
  /* Be a good citizen and clean up... */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dvector(V,0,Nmx);
  free_dmatrix(E,0,lmax,0,nmaxmax);
  free_dmatrix(F,0,lmax,0,nmaxmax);
  free_d3tensor(Psi,0,lmax,0,nmaxmax,0,Nmx);
  free_dvector(Rho,0,Nmx);
  free_dvector(Rhonew,0,Nmx);
  free_dvector(phi,0,Nmx);
  free_dvector(Vxc,0,Nmx);
  free_dvector(Depsxc,0,Nmx);
  free_dvector(integrand,0,Nmx);
  free_ivector(nmax,0,lmax);
  
}
