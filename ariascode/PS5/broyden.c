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
  double **Rho,**g,Z;
  int lmax,*nmax,nmaxmax;
  double **G,*u;
  double **F,**LU,*b;
  double Etot,sgn,tmp;
  int *indx, itp, itpp,i,j;

  LU=dmatrix(1,Itmx+2,1,Itmx+2); /* Space for LU decomposition of matrix */
  b=dvector(1,Itmx+2); /* Space for answer to G*b=u */
  indx=ivector(1,Itmx+2); /* Space for integer vector needed by ludcmp */


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
  Rho=dmatrix(1,Itmx+1,0,Nmx);
  g=dmatrix(1,Itmx+1,0,Nmx);
  G=dmatrix(1,Itmx+2,1,Itmx+2);
  u=dvector(1,Itmx+2);

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
    Rho[1][k]=0.;
  
  for (it=1; it<=Itmx; it++) {
    Etot=getg(g[it],Rho[it],Z,lmax,nmax,nmaxmax,F,r,dr,N);
    printf("%4d  %20.12f  %20.12f\n",
           it,fabs(Etot-(-74.473076803203738)),log(fabs(Etot-(-74.473076803203738))));
    
    for (k=0; k<=N; k++)
      Rho[it+1][k]=Rho[it][k]+alpha*g[it][k];
    getg(g[it+1],Rho[it+1],Z,lmax,nmax,nmaxmax,F,r,dr,N);
    for (itp=1;itp<=it+1;itp++)
      for (itpp=1;itpp<=it+1; itpp++){
    	G[itp][itpp]=0.;
	for(k=0;k<=N;k++)
	  G[itp][itpp]+=g[itp][k]*g[itpp][k];
      }
    for(k=1;k<=it+1;k++)
      G[it+2][k] = 1.;
    for(k=1;k<=it+1;k++)
      G[k][it+2] = -0.5;
    G[it+2][it+2]=0.;
    for(k=1;k<=it+1;k++)
      u[k]=0;
    u[it+2]=1;
    
    for (i=1; i<=it+2; i++) /* Copy G into LU because ludcmp destroys its input */
      for (j=1; j<=it+2; j++)
    	LU[i][j]=G[i][j];
    for (i=1; i<=it+2; i++) /* Copy u into b because lubksb destroys its input */
      b[i]=u[i];
 
    ludcmpp480(LU,it+2,indx,&sgn); /* Get LU decomposition */
    lubksbp480(LU,it+2,indx,b); /* Use LU decomposition to solve equations */
    
    //        printf("\nTesting G*b=u?:\n");
    //	for (i=1; i<=it+2; i++) {
    //	  tmp=0.;
    //	  for (j=1; j<=it+2; j++)
    //	    tmp+=G[i][j]*b[j];
    //	  printf("%20.12f  =%20.12f  ?\n",tmp,u[i]);
    //	}
    //   printf("\n");
   
   //    for (k=0;k<=N;k++)
   // for(itp=it+1;itp>=1;itp--)
   //   Rho[it+1][k]+=Rho[itp][k]*b[itp];
   
   for (k=0; k<=N; k++)
     Rho[it+1][k]*=b[it+1];

   for (k=0; k<=N; k++)
     for (itp=1;itp<=it;itp++)
       Rho[it+1][k]+=b[itp]*Rho[itp][k];

   for (k=0;k<=N; k++)
     Rho[it+1][k]=fabs(Rho[it+1][k]);

 
  }
  
  /* Be a good citizen and clean up... */
  free_dvector(r,0,Nmx);
  free_dvector(dr,0,Nmx);
  free_dmatrix(G,1,Itmx+2,1,Itmx+2);
  free_dvector(u,1,Itmx+2);
  free_dmatrix(LU,1,Itmx+2,1,Itmx+2);
  free_dvector(b,1,Itmx+2);
  free_ivector(indx,1,Itmx+2);

  free_dmatrix(F,0,lmax,0,nmaxmax);
  free_dmatrix(Rho,1,Itmx+1,0,Nmx);
  free_dmatrix(g,1,Itmx+1,0,Nmx);
}
