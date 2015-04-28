#include <stdio.h>
#include "p480.h"
#include <math.h>
#include "nrutil.h"

#define JMAX 80
#define MAXIT 120
#define UNUSED (-1.11e30)
#define NR_END 1
#define FREE_ARG char*


#define NRANSI
#define TINY 1.0e-20;

void ludcmpp480(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY
#undef NRANSI


void lubksbp480(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}




double ***d3tensor(long nrl,long nrh,long ncl,long nch, long ndl, long ndh){
  /* allocate a double 3tensor with range t[nrl..nrh][ncl..ncl][ndl..ndh] */
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
        double ***t;

        /* allocate pointers to pointers to rows */
        t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
        if (!t) nrerror("allocation failure 1 in d3tensor()");
        t += NR_END;
        t -= nrl;

        /* allocate pointers to rows and set pointers to them */
        t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
        if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
        t[nrl] += NR_END;
        t[nrl] -= ncl;

        /* allocate rows and set pointers to them */
        t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
        if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
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

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh)
/* free a double d3tensor allocated by f3tensor() */
{
        free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
        free((FREE_ARG) (t[nrl]+ncl-NR_END));
        free((FREE_ARG) (t+nrl-NR_END));
}

  

void interpolate(double f[], int N){
  int k;
  f[1]=0.0625*(9*(f[0]+f[2])-f[4]);
  f[N-1]=0.0625*(-f[N-4]+9*(f[N-2]+f[N]));
  for (k=3;k < N-1; k+=2)
    f[k]=0.0625*(-f[k-3]+9.0*(f[k-1]+f[k+1])-f[k+3]);
}

double zriddrp480(double (*func)(double,int,double[],double[],double[],int),
                  double x1, double x2, double xacc,
                  int match, double V[], double r[], double dr[], int N){
  
  int j;
  double ans, fh,fl,fm,fnew,s,xh,xl,xm,xnew;

  fl=(*func)(x1,match,V,r,dr,N);
  fh=(*func)(x2,match,V,r,dr,N);
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl=x1;
    xh=x2;
    ans=UNUSED;
    for (j=1;j<=MAXIT;j++) {
      xm=0.5*(xl+xh);
      fm=(*func)(xm,match,V,r,dr,N);
      s=sqrt(fm*fm-fl*fh);
      if (s == 0.0) {
	//printf("(%d iterations in zriddrp480)\n",j);
	return ans;
      }
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
      if (fabs(xl-xh) <= xacc) return ans;
      ans=xnew;
      fnew=(*func)(ans,match,V,r,dr,N);
      if (fnew == 0.0) return ans;
      if (SIGN(fm,fnew) != fm) {
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
      if (fabs(xl-xh) <=xacc) {
	//	//printf("(%d iterations in zriddrp480, func=%e)\n",j,fl);
	return ans;
      }
    }
    //nrerror("zriddrp480 exceed maximum iterations");
  }
  else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    //nrerror("root must be bracketed in zriddrp480.");
  }
  return 0.0;
}

      
double rtbisp480(double (*func)(double,int,double [],double [],double[],int),
                 double x1, double x2, double xacc,
                 int match, double V[], double r[], double dr[], int N){

  void nrerror(char error_text[]);
  int j;
  double dx,f,fmid,xmid,rtb;

  f=(*func)(x1,match,V,r,dr,N);
  fmid=(*func)(x2,match,V,r,dr,N);

  if (f*fmid > 0.0) //nrerror("root must be bracketed for bisection in rtbisp480");
  if (f == 0.0) {
    //printf("(0 iterations in rtbisp480, func=%e)\n",f);
    return x1;
  }
  if (fmid == 0.0) {
    //printf("(0 iterations in rtbisp480, func=%e)\n",fmid);
    return x2;
  }
  
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); 
  /* Bisection loop */
  for (j=1;j<=JMAX;j++) {
    fmid=(*func)(xmid=rtb+(dx *= 0.5),match,V,r,dr,N);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0){
      //printf( "  %d iterations in rtbisp480\n",j);
      return rtb;
    }

  }
  // nrerror("Too many bisections in rtbis480");
  return 0.0;   /* Never get here. */
}

void rk4p480(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double [],int, double [], double [], double [], int), int k, int dk, double f[], double r[], double dr[], int N)
{
  int i;
  double xh,hh,h6,kh,dym[3],dyt[3],yt[3];
  //double *dym, *dyt, *yt;
  /*
  dym=dvector(1,n);
  dyt=dvector(1,n);
  yt=dvector(1,n);
  */
  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  kh=k+dk*0.5;
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
  (*derivs)(xh,yt,dyt,kh,f,r,dr,N);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(xh,yt,dym,kh,f,r,dr,N);
  for (i=1;i<=n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt,k+dk,f,r,dr,N);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  /*
  free_dvector(yt,1,n);
  free_dvector(dyt,1,n);
  free_dvector(dym,1,n);
  */
}

double simpint(double f[], double r[], double dr[], int N) {
        
  /*Declare working variables*/
  double Trap1,Trap2;
  double Total=0;
  double Total2=0;
  double Simp;
  int i,j;

  /* Calculate value of each interval and add them up for spacing=h */
  for (i=1;i<=N; i++) {
    Trap1=(f[i-1]*dr[i-1]+f[i]*dr[i])/2/N;
    Total=Total+Trap1;
  }

  /* Calculate value of each interval and add them up for spacing=2h */
  for (j=2; j<=N; j=j+2) {
    Trap2=(f[j-2]*dr[j-2] + f[j]*dr[j])/N;
    Total2 = Total2 + Trap2;
  }
  
  /* Richardson Extrapolation */
  Simp=(((double)4)/((double)3))*Total-(((double)1)/((double)3))*Total2; 
  return Simp;
}

double trapint(double f[], double r[], double dr[], int N) {
  
  /* Declare working variables */ 
  double Interval;
  double Trap=0;
  int i;

  /* Calculate value of each interval and add them up for spacing=h */
  for (i=1; i<=N; i++) {
    Interval=(f[i-1]*dr[i-1]+f[i]*dr[i])/(2*N);
    Trap=Trap+Interval;
  }
  return Trap;
}
