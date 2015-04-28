#include <stdio.h>
#include <math.h>
#include "physics.h"
#include "p480.h"
#include "nrutil.h"

main()
{
  /* Grid information */
  double *r,*dr;
  int N;

  /* Physics vectors */
  double *V;

  /* Working variables */
  int k,Nstates;
  double x,Elower;
  double *E;
  
  /* Value of pi */
  double pi;
  pi=4.*atan(1.);

  /* Loop over different input N and numbers of states */
  while (1) {
    printf("N="); scanf("%d",&N);
    printf("Nstates="); scanf("%d",&Nstates);

    /* Allocate NR vectors */
    r=dvector(0,N);
    dr=dvector(0,N);
    V=dvector(0,N);

    E=dvector(0,Nstates-1);

    /* Construct change of variables information and
       integrand function for computing V_el-nuc */
    for (k=0; k<=N; k++) {
      x=((double) k)/((double) N);
      r[k]=x*x/(1.-x);
      dr[k]=1./(1.-x)/(1.-x)-1.;
      V[k]=-1./r[k];
    }
    /* Special conditions for handling end-points of integration */
    V[0]=0.;
    dr[N]=0.;


    getEs(E,Nstates-1,Elower=-1.,V,r,dr,N);

    for (k=0; k<Nstates; k++)
      printf("State with %2d nodes, E=%16.12f\n",k,E[k]);

    /* Clean up vectors before next set is allocated! */
    free_dvector(r,0,N);
    free_dvector(dr,0,N);
    free_dvector(V,0,N);
    free_dvector(E,0,Nstates-1);
  }

}


