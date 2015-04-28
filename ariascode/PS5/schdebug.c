#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

int main()
{
  /* Grid information */
  double *r,*dr;
  int N;

  /* Physics vectors */
  double *V,*Psianal,*Psiout;

  /* Working variables */
  int k;
  double x;
  int k1,k2,nodes;
  double Psi,Psip,E;

  /* Value of pi */
  double pi;
  pi=4.*atan(1.);

  /* Loop over different input energies for fixed N=100 */
  N=100;
  while (1) {
    printf("E=");
    scanf("%lf",&E);

    /* Allocate NR vectors */
    r=dvector(0,N);
    dr=dvector(0,N);
    V=dvector(0,N);
    Psianal=dvector(0,N);
    Psiout=dvector(0,N);

    /* Construct change of variables information and
       integrand function for computing V_el-nuc */
    for (k=0; k<=N; k++) {
      x=((double) k)/((double) N);
      r[k]=x*x/(1.-x);
      dr[k]=1./(1.-x)/(1.-x)-1.;
      V[k]=-1./r[k];
      Psianal[k]=r[k]*exp(-r[k]);
    }
    /* Special conditions for handling end-points of integration */
    V[0]=0.;
    dr[N]=0.;

    /* Outward integration test... */
    Psi=0.;
    Psip=1.;
    /* Because Psi and Psip contain both input and output values, they
       must be passed as pointers (hence the '&') */
    /* Psiout, as an array, is already a pointer, so no '&' is needed. */
    /* We could just use "0,N," instead of "k1=0,k2=N,", but this notation
       makes the code easier to read, i.e., debug. */
    nodes=schint(&Psi,&Psip,Psiout,k1=0,k2=N,V,E,r,dr,N);

    /* Only compare on even points, as these are the only places where
       RK generates values */

    printf("=============================================================\n");
    printf("Outward integration test:\n");
    printf("-------------------------\n");
    printf("   %d nodes\n",nodes);
    printf("   Return value of Psi: %16.12f\n\n",Psi);

    printf("   k           r             Psianal          Psiout       Psiout/Psianal\n");
    for (k=0; k<=N; k+=2) {
      printf("%5d %16.12f %16.12f %16.12f %16.12f\n",
             k,r[k],Psianal[k],Psiout[k],Psiout[k]/Psianal[k]);
    }

    /* Inward integration test... */
    Psi=0.;
    Psip=-1.;
    /* Because Psi and Psip contain both input and output values, they
       must be passed as pointers (hence the '&') */
    /* Psiout, as an array, is already a pointer, so no '&' is needed. */
    /* We could just use "0,N," instead of "k1=0,k2=N,", but this notation
       makes the code easier to read, i.e., debug. */

    nodes=schint(&Psi,&Psip,Psiout,k1=N,k2=0,V,E,r,dr,N);

    /* Only compare on even points, as these are the only places where
       RK generates values */
    
    printf("=============================================================\n");
    printf("Inward integration test: %d nodes\n",nodes);
    printf("-------------------------\n");
    printf("   Return value of Psi: %16.12e\n\n",Psi);

    printf("   k           r             Psianal          Psiout       Psiout/Psianal\n");
    for (k=0; k<=N; k+=2) {
      printf("%5d %16.12f %16.12f %16.12e %16.12e\n",
             k,r[k],Psianal[k],Psiout[k],Psiout[k]/Psianal[k]);
    }

    /* Clean up vectors before next set is allocated! */
    free_dvector(r,0,N);
    free_dvector(dr,0,N);
    free_dvector(V,0,N);
    free_dvector(Psianal,0,N);
    free_dvector(Psiout,0,N);
    
  }
  
}

