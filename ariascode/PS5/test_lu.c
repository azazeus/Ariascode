#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "p480.h"
#include "physics.h"

#define N 4
main()
{
  double **G,**LU,*b,*u;

  int *indx;
  double sgn,tmp;

  int i,j;

  /* Declare space */
  G=dmatrix(1,N,1,N); /* Matrix */
  LU=dmatrix(1,N,1,N); /* Space for LU decomposition of matrix */
  b=dvector(1,N); /* Space for answer to G*b=u */
  u=dvector(1,N); /* Space for right-hand side, u */
  indx=ivector(1,N); /* Space for integer vector needed by ludcmp */

  /* Make matrix Gij=1/(i+j) */
  for (i=1; i<=N; i++)
    for (j=1; j<=N; j++)
      G[i][j]=1./( (double) i+j );

  /* Make right-hand side u=(1 1 1 ...) for equation G*x=u */
  for (i=1; i<=N; i++) /* u=0 for all but the last component */
    if (i<N)
      u[i]=0.;
    else
      u[i]=1.;


  /* Print out the problem */
  printf("\nSolving...\n\n");
  for (i=1; i<=N; i++) {
    printf("| ");
    for (j=1; j<=N; j++)
      printf("%9.2e ",G[i][j]);
    printf("|   |b[%d]|  =  | %9.2e |\n",i,u[i]);
  }

  /* Solve G*b=u */
  for (i=1; i<=N; i++) /* Copy G into LU because ludcmp destroys its input */
    for (j=1; j<=N; j++)
      LU[i][j]=G[i][j];
  for (i=1; i<=N; i++) /* Copy u into b because lubksb destroys its input */
    b[i]=u[i];

  /* Actual solution of equations */
  ludcmpp480(LU,N,indx,&sgn); /* Get LU decomposition */
  lubksbp480(LU,N,indx,b); /* Use LU decomposition to solve equations */

  /* Does it work? */
  printf("\nTesting G*b=u?:\n");
  for (i=1; i<=N; i++) {
    tmp=0.;
    for (j=1; j<=N; j++)
      tmp+=G[i][j]*b[j];
    printf("%20.12f  =%20.12f  ?\n",tmp,u[i]);
  }
  printf("\n");

  /* Be good, and clean up! */
  free_dmatrix(G,1,N,1,N); /* Matrix */
  free_dmatrix(LU,1,N,1,N); /* Space for LU decomposition of matrix */
  free_dvector(b,1,N); /* Space for answer to G*b=u */
  free_dvector(u,1,N); /* Space for right-hand side, u */
  free_ivector(indx,1,N); /* Space for integer vector needed by ludcmp */
}
