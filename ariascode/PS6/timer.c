#include <stdio.h>
#include <sys/time.h>
#include "nrutil.h"
#include <math.h>
#include "perform.h"

#define Nb 32

double mult(int n,double **M1,double **M2,double **M3,
	    double *v1,double *v2,double *v3)
{
  int i;
  double x=1.;
  
  for (i=0; i<n*n*n; i++)
    x=x*1.000001;

  return x;
}

double divide(int n,double **M1,double **M2,double **M3,
            double *v1,double *v2,double *v3)
{
  int i;
  double x=1.;

  for (i=0; i<n*n*n; i++)
    x=x/1.000001;

  return x;
}

double multsub(int n,double **M1,double **M2,double **M3,
            double *v1,double *v2,double *v3)
{
  int i;
  double x=1.;
  double y=1.000001;
  for (i=0; i<n*n*n; i++)
    //x=x*1.000001;
    x=sub(x);

  return x;
}

double multif(int n,double **M1,double **M2,double **M3,
            double *v1,double *v2,double *v3)
{
  int i;
  double x=1.;

  for (i=0; i<n*n*n; i++){
    if(x<1.000001){
      x=x*1.000001;
    }else{x=x*0.99;}
  }
  return x;
}

double multal(int n,double **M1,double **M2,double **M3,
            double *v1,double *v2,double *v3)
{
  int i;
  double x=1., *p;
  for (i=0; i<n*n*n; i++){
    p=dvector(0,1);
    x=x*1.000001;
    free_dvector(p,0,1);
  }
  return x;
}

double multadd(int n,double **M1,double **M2,double **M3,
            double *v1,double *v2,double *v3)
{
  int i;
  double x=1.;

  for (i=0; i<n*n*n; i++)
    x+=x*.000001;

  return x;
}

double cmat1(int n,double **M1,double **M2,double **M3,
	       double *v1,double *v2,double *v3)
{
  int i,j,k;
  double x=1.;

  for (k=0; k<n; k++)
    for (i=0; i<n;i++)
      for (j=0;j<n; j++)
	M3[i][j]=M1[i][j];

  return x;
}

double cmat2(int n,double **M1,double **M2,double **M3,
	     double *v1,double *v2,double *v3)
{
  int i,j,k;
  double x=1.;

  for (k=0; k<n; k++)
    for (j=0; j<n;j++)
      for (i=0;i<n; i++)
        M3[i][j]=M1[i][j];

  return x;
}

double vecmat1(int n,double **M1,double **M2,double **M3,
             double *v1,double *v2,double *v3)
{
  int i,j,k;
  double x=1.;

  for (k=0; k<n; k++)
    for (i=0; i<n;i++)
      for (j=0;j<n; j++)
        v2[i]=v1[j]*M1[j][i];

  return x;
}

double vecmat2(int n,double **M1,double **M2,double **M3,
	       double *v1,double *v2,double *v3)
{
  int i,j,k;
  double x=1.;

  for (k=0; k<n; k++)
    for (j=0; j<n;j++)
      for (i=0;i<n; i++)
        v2[i]=v1[j]*M1[j][i];

  return x;
}

double matmat1(int n,double **M1,double **M2,double **M3,
               double *v1,double *v2,double *v3)
{
  int i,j,k;
  double x=1.;

  for (i=0; i<n; i++)
    for (j=0; j<n;j++)
      for (k=0;k<n; k++)
        M3[i][j]+=M1[i][k]*M2[k][j];

  return x;
}

double matmat2(int n,double **M1,double **M2,double **M3,
               double *v1,double *v2,double *v3)
{
  int i,j,k,I,J,K,p,q;
  double x=1.;
  double b1[Nb][Nb], b2[Nb][Nb],b3[Nb][Nb];

  for (I=0; I<n/Nb; I++){
    for (J=0; J<n/Nb;J++){
      for (p=0;p<Nb;p++)
	for(q=0;q<Nb;q++)
	  b3[p][q]=0.;
      for(K=0;K<n/Nb;K++){
	/* Copy IK block of M1 into b1, an Nb x Nb matrix */
	for (i=0; i<Nb; i++) 
	  for (k=0; k<Nb; k++)
	    {b1[i][k]=M1[I*Nb+i][K*Nb+k];
	    }

	/* Copy KJ block of M2 into b2, an Nb x Nb matrix */
	for (i=0; i<Nb; i++)
          for (k=0; k<Nb; k++)
            {b2[i][k]=M2[K*Nb+i][J*Nb+k];
	    }
	/* Do b1*b2, adding result into b3 */
	for (i=0; i<Nb; i++)
	  for (j=0; j<Nb; j++)
	    for (k=0; k<Nb; k++)
	      {b3[i][j]+=b1[i][k]*b2[k][j];
	      }
	
	      }
      /* Copy IJ block of b3 into M3, an Nb x Nb matrix */
      for (i=0; i<Nb; i++)
	for (k=0; k<Nb; k++)
	  {M3[I*Nb+i][J*Nb+k]=b3[i][k];
	  }
    }
  }
  return x;
}
  

double timer(char *name,double flop,
             double (*func)(int,double **,double **,double **,
                            double *,double *,double *),
             int n)
{
  /* Timer variables */
  struct timeval time1,time2;
  double ctime,mflop;
  double **M1,*v1;
  double **M2,*v2;
  double **M3,*v3;
  int i,j;

  M1=dmatrix(0,n,0,n);
  M2=dmatrix(0,n,0,n);
  M3=dmatrix(0,n,0,n);

  v1=dvector(0,n);
  v2=dvector(0,n);
  v3=dvector(0,n);

  /* Fill in values for input vectors */
  for (i=0; i<n; i++) {
    v1[i]=3.14159;
    v2[i]=0.314159;
  }

  /* Fill in values for input matrices */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      M1[i][j]=1./(1+i+j);
      M2[i][j]=1./(2+i+j);
      M3[i][j]=0.;
    }

  /* Record time of function call */
  gettimeofday(&time1,NULL);

  /* Call function */
  func(n,M1,M2,M3,v1,v2,v3);

  /* Record time of function return */
  gettimeofday(&time2,NULL);

  /* Coupute time in seconds */
  ctime=
    (time2.tv_sec-time1.tv_sec) /* Integer seconds */
    +(time2.tv_usec-time1.tv_usec)/1.e6; /* Integer giving microseconds part */

  /* Convert from FLOP to MFLOP */
  mflop=flop/1e6;

  /* Output result */
  fprintf(stderr,"%20s: %4.0lf MFLOPS  = %.1lf MFLOP / %.3f sec\n",
          name,mflop/ctime,mflop,ctime);

  /* Free memory before return */
  free_dvector(v1,0,n);
  free_dvector(v2,0,n);
  free_dvector(v3,0,n);

  free_dmatrix(M1,0,n,0,n);
  free_dmatrix(M2,0,n,0,n);
  free_dmatrix(M3,0,n,0,n);

  return ctime;
}

main()
{
  /* Place to store run times of various routines */
  double tm,td,ts,ti,tal,tma,tac1,tac2,vm1,vm2;

  /* Variables for verification of matmat3 */
  int i,j;
  double **M1,**M2,**M3,**M4;
  double err=0.;

  /* Problem size (n x n) */
  int n;

  printf("Problem size: ");
  scanf("%d",&n);

  /* Calls to timers which don't access much memory */
  tm=timer("mult",(double) n*n*n,mult,n);
  td=timer("divide",(double) n*n*n,divide,n);
  ts=timer("sub ",(double) n*n*n,multsub,n); 
  ti=timer("if ",(double) n*n*n,multif,n); 
  tal=timer("alloc",(double) n*n*n,multal,n); 
  tma = timer("multadd",(double) 2*n*n*n,multadd,n); 

  /* Time strides through memory */
  tac1=timer("cmat1",(double) n*n*n,cmat1,n);
  tac2=timer("cmat2",(double) n*n*n,cmat2,n); 

  /* Print out time ratios */
  printf("Notes:\n"); 
  printf("Division      costs %4.1f mults\n",td/tm); 
  printf("Subroutine    costs %4.1f mults\n",ts/tm-1);
  printf("If            costs %4.1f mults\n",ti/tm-1); 
  printf("malloc        costs %4.1f mults\n",tal/tm-1);
  printf("multadd       costs %4.1f mults\n",tma/tm);

  printf("Memory access costs %4.1f mults\n",tac1/tm-1);
  printf("Cache miss    costs %4.1f accesses\n",tac2/tac1);
  printf("\n");

  /* Matrix-vector operations */
  timer("vecmat1",(double) 2*n*n*n,vecmat1,n);
  timer("vecmat2",(double) 2*n*n*n,vecmat2,n);
  timer("matmat1",(double) 2*n*n*n,matmat1,n); 
  timer("matmat2",(double) 2*n*n*n,matmat2,n); 

  /* Verify matmat2 */
    M1=dmatrix(0,n,0,n); 
    M2=dmatrix(0,n,0,n); 
    M3=dmatrix(0,n,0,n); 
    M4=dmatrix(0,n,0,n); 
  
  /* Fill in values for input matrices */
    for (i=0; i<n; i++) 
      for (j=0; j<n; j++) { 
	M1[i][j]=1./(1+i+j); 
	  M2[i][j]=1./(2+i+j); 
      } 

  /* Note, the matmat routines don't use their last three arguments anyway */
    matmat1(n,M1,M2,M3,M1[0],M2[0],M3[0]); 
    matmat2(n,M1,M2,M4,M1[0],M2[0],M3[0]); 

    for (i=0; i<n; i++) 
      for (j=0; j<n; j++) 
	err+=pow(M4[i][j]-M3[i][j],2);
    printf("\nIMPORTANT: Test of matmat2; should be near zero... %e\n\n",err); 

    /* Free up space */
    free_dmatrix(M1,0,n,0,n); 
    free_dmatrix(M2,0,n,0,n);
    free_dmatrix(M3,0,n,0,n);
    free_dmatrix(M4,0,n,0,n);
}

