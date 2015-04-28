#include <math.h>
#include <stdio.h>
double excPZ(double rs);

double excpPZ(double rs);

main()
{
  double rs,drs;

  rs=3.; /* Test rs > 1 "branch" */
  printf("Testing derivative for rs=%f ...\n",rs);
  for (drs=1; drs>1e-6; drs/=10)
    printf("    %20.16f %20.16f\n",drs,
           (excPZ(rs+drs)-excPZ(rs))/drs /* <- Finite difference slope */
           /                             /* Ratio should approach 1! */
           excpPZ(rs)                    /* <- Coded value of derivative */
           );

  printf("\n");

  rs=0.3; /* Test rs < 1 "branch" */
  printf("Testing derivative for rs=%f ...\n",rs);
  for (drs=1; drs>1e-6; drs/=10)
    printf("    %20.16f %20.16f\n",drs,
           (excPZ(rs+drs)-excPZ(rs))/drs /* <- Finite difference slope */
           /                             /* Ratio should approach 1! */
           excpPZ(rs)                    /* <- Coded value of derivative */
           );
}

