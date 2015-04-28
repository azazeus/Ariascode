#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include "nrutil.h"
#define NR_END 1
#define FREE_ARG char*


double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
    /* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)))\
    ;
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
     /* free a double d3tensor allocated by d3tensor() */
{
free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
free((FREE_ARG) (t[nrl]+ncl-NR_END));
free((FREE_ARG) (t+nrl-NR_END));
}



void rk4p480(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double [],int, double [], double [], double [], int), int k, int dk, double f[], double r[], double dr[], int N )
{
  
  int i;
  double xh,hh,h6,hk,kh,*dym,*dyt,*yt;
  dym=dvector(1,n);
  dyt=dvector(1,n);
  yt=dvector(1,n);
  hk=dk*0.5;
  h6=h/6.0;
  kh=k+hk;
  
  for (i=1;i<=n;i++) yt[i]=y[i]+(h/2.)*dydx[i]; // First step.
  (*derivs)(kh,yt,dyt,kh,f,r,dr,N); //Second step.
  for (i=1;i<=n;i++) yt[i]=y[i]+(h/2.)*dyt[i];
  (*derivs)(kh,yt,dym,kh,f,r,dr,N); //Third step.
  for (i=1;i<=n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(k+h,yt,dyt,(k+dk),f,r,dr,N); //Fourth step.
  for (i=1;i<=n;i++) //Accumulate increments with proper weights. 
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  free_dvector(yt,1,n);
  free_dvector(dyt,1,n);
  free_dvector(dym,1,n);
}

double trapint(double f[], double r[], double dr[], int N)
{
  int l=0;
  double sum = 0.0;
  for(l=1; l<=N-1; l++){
    sum = sum + f[l]*dr[l];
  }
  sum = sum + 0.5*(f[0]*dr[0] + f[N]*dr[N]);
  return (double)(sum/N);
}

double simpint(double f[], double r[], double dr[], int N)
{
  int l=0;
  double sum = 0.0, s=4.;
  for(l=1; l<=N-1; l++){
    sum = sum + s*f[l]*dr[l];
    s=6.-s;
  }
  sum = sum + f[0]*dr[0] + f[N]*dr[N];
  return (double)(sum/(3.0* (double)N));
}

#define JMAX 80 //Maximum allowed number of bisections.

double rtbisp480(double (*func)(double,int,double [], double [], double [], int), double x1, double x2, double xacc, int match, double V[], double r[], double dr[], int N)
{
  void nrerror(char error_text[]);
  int j;
  double dx,f,fmid,xmid,rtb;
  f=(*func)(x1,match,V,r,dr,N);
  fmid=(*func)(x2,match,V,r,dr,N);
  if (f*fmid > 0.0) nrerror("Root must be bracketed for bisection in rtbisp480");
  if(f==0){
    printf("(0 iterations in rtbisp480, func=%e)\n",f);
    return x1;
    exit(0);
  }
  if(fmid==0){
    printf("(0 iterations in rtbisp480, func=%e)\n",fmid);
    return x2;
    exit(0);
  }
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); 
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5),match,V,r,dr,N); //bisection loop
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) {
      printf("   (%d iterations in rtbisp480, func =%e)\n", j, fmid);
      return rtb;
    }
  }
  nrerror("Too many bisections in rtbis");
  return 0.0; //Never get here.
  
}

#define MAXIT 60
#define UNUSED (-1.11e30)

double zriddrp480(double (*func)(double,int,double [],double [], double [],int), double x1, double x2, double xacc, int match, double V[], double r[], double dr[], int N)
{
  int j;
  double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
  fl=(*func)(x1,match,V,r,dr,N);
  fh=(*func)(x2,match,V,r,dr,N);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=UNUSED; //Any highly unlikely value, to simplify logic below.
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm=(*func)(xm,match,V,r,dr,N); /*First of two function evaluations per iteration. */
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0){ 
	printf("  (%d iterations in zriddrp480, func = %e)\n",j,fl);
	return ans;
      }
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); //Updating formula.
      if (fabs(xl-xh) <= xacc) {
	printf("  (%d iterations in zriddrp480, func = %e)\n",j,fl);
	return ans;
      }
      ans=xnew;
      fnew=(*func)(ans,match,V,r,dr,N); /*Second of two function evaluations per iteration. */
      if (fnew == 0.0) {
	printf("  (%d iterations in zriddrp480, func = %e)\n",j,fnew);
	return ans;
      }
      if (SIGN(fm,fnew) != fm) 
	{// Bookkeeping to keep the root bracketed on next iteration. 
	  xl=xm;
	  fl=fm;
	  xh=ans;
	  fh=fnew;
	} else if (SIGN(fl,fnew) != fl) {
	  xh=ans;
	  fh=fnew;
	} else if (SIGN(fh,fnew) != fh) {
	  xl=ans;
	  fl=fnew;
	} else nrerror("never get here.");
      if (fabs(xh-xl) <= xacc) {
	printf("  (%d iterations in zriddrp480, func = %e)\n",j,fl);
	return ans;
      }
    }
    nrerror("zriddr exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    nrerror("root must be bracketed in zriddr.");
  }
  return 0.0; //Never get here.
}

void interpolate(double f[], int N)
{
  int k;
  for(k=3; k<=N-3; k+=2){
    f[k] = 0.0625*(-f[k-3] + 9*(f[k-1]+f[k+1]) - f[k+3]);
  }
  f[1]=0.0625*(9*(f[0]+f[2]) - f[4]);
  f[N-1]=0.0625*(-f[N-4] + 9*(f[N-2]+f[N]));
}



