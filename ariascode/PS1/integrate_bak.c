/*
  Main program for PS#1
  Tests Simpson and Trapezoid integration routines
*/


/* Standard C header for using I/O and math libraries */
#include <stdio.h>
#include <math.h>

double trapint(double f[], double r[], double dr[], int N)
{
  int l=0;
  double sum = 0.0;
  for(l=1; l<=N-1; l++){
    sum = sum + f[l];
  }
  sum = sum + 0.5*(f[0] + f[N]);
  return (sum/(double)N);
}
  
double simpint(double f[], double r[], double dr[], int N)
{
  int l=0;
  int k=2;
  double sum = 0.5*(f[0]+f[N]);
  for(l=1; l<=N-1; l++){
	  k=6-k;
	  sum=(((double)l+1.0)/((double)l+2.0))*sum+((double)k*f[l]/((double)l+2.0));
    //sum = sum + 2.0*(l%2 + 1.0)*f[l];
  }
  //sum = sum + f[0] + f[N];
  return ((sum*((double)N+1.0)/((double)N))/3.0);
}
  
  /* Maximum size of N for integrations */
#define Nmx 200000
  main()
    {
      /* Declare space for arrays, remember that there are N+1 points! */
      /* Double precision for accurate numerical work... */
      double r[Nmx+1],dr[Nmx+1],f[Nmx+1];
      
      /* Declare working variables */
      int i,N;
      double u;
      
      /* Trick to get all digits of PI to the full machine precision */
      double pi;
      pi=4.*atan(1.);
      
      
      /* Loop over sizes 2, 20, 200, 2000, ..., Nmx */
      for (N=2; N<=Nmx; N*=10) {
	/* For each size, fill in arrays */
	for (i=0; i<=N; i++) {
	  /* Let u be equally spaced between 0 and 1, -inclusive- */
	  u=((double) i)/((double) N);
	  
	  /* Define change of variables */
	  r[i]=u;
	  dr[i]=1.;
	  
	  /* Specify function */
	  f[i]=1/(1+r[i]*r[i])/(1+r[i]*r[i])/(0.25+0.125*pi);
	}
	
	/* Special cases at end points (reasons for this will become clear...) */
	f[0]=f[0];
	r[0]=r[0];
	dr[0]=dr[0];
	
	f[N]=f[N];
	r[N]=r[N];
	dr[N]=dr[N];
	
	/* Call integration routines and output results... */
	printf("Trap, Simp: %7d %.14E,%.14E\n",N,1.0-trapint(f,r,dr,N),1.0-simpint(f,r,dr,N));
      }  
      
    }
    
