#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"
#include "physics.h"

void derivs_Schrodinger(double, double *, double *, int, double *, double *, double *, int);

int schint(double *,double *,double *,int,int,double *,double,double *,double *,int);

double func_Schrodinger(double, int, double *, double *, double *, int);

double func_SchrodingerNodes(double, int, double *, double *, double *, int);

void getES(double *, int, double, double *, double *, double *, int);

void getallEs(double **, int, int *, double, double *, double *, double *, int);

void getPsi(double *, double, int, double *, double *, double *,int);

void getallPsis(double ***, double **,int, int *,double *,double *,double *, int);

void getRho(double *, double ***, double **, int, int *, int);

void getphi(double *, double *, double *, double *, int);

void derivs_poisson(double, double *, double *, int, double *, double *, double *, int);

double excPZ(double);

double excpPZ(double);

void getVxc(double *, double *, double *, double *, int);

void getDepsxc(double *, double *, double *, double *, int);

double exc(double);

double excp(double);

double getg(double *, double *,double, int, int *, int, double **, double *, double *, int);

void rk4p480sch(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double [],int, double [], double [], double [], int), int k, int dk, double f[], double r[], double dr[], int N );


/* ============ end of function declarations================================ */

void rk4p480sch(double y[], double dydx[], int n, double x, double h, double yout[], void (*derivs)(double, double [], double [],int, double [], double [], double [], int), int k, int dk, double f[], double r[], double dr[], int N )
{
  
  int i;
  double xh,hh,h6,hk,kh,dym[3],dyt[3],yt[3];
  /*dym=dvector(1,n);
  dyt=dvector(1,n);
  yt=dvector(1,n);
  hk=dk*0.5;*/
  h6=h/6.0;
  kh=k+hk;
  
  for (i=1;i<=n;i++) yt[i]=y[i]+(h/2.)*dydx[i]; // First step.
  derivs_Schrodinger(kh,yt,dyt,kh,f,r,dr,N); //Second step.
  for (i=1;i<=n;i++) yt[i]=y[i]+(h/2.)*dyt[i];
  derivs_Schrodinger(kh,yt,dym,kh,f,r,dr,N); //Third step.
  for (i=1;i<=n;i++) {
    yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(k+h,yt,dyt,(k+dk),f,r,dr,N); //Fourth step.
  for (i=1;i<=n;i++) //Accumulate increments with proper weights. 
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  /*  free_dvector(yt,1,n);
  free_dvector(dyt,1,n);
  free_dvector(dym,1,n);*/
}

void derivs_Schrodinger(double x, double y[], double dydx[], int k, double T[], double r[], double dr[], int N){
  dydx[1] = y[2]*dr[k];
  dydx[2] = -2.*T[k]*y[1]*dr[k];
}


#define SMALL (1.e-300)
#define SQRTSMALL (1.e-150)
#define BIG (1.e300)
#define SQRTBIG (1.e150)

int schint(double *Psi, double *Psip, double Psiout[], int k1, int k2, double V[], double E, double r[], double dr[], int N){
  
  double *y, *dydx, *T ;
  int i, k, g;
  int nodes = 0;
  double pi;
  pi = 4.*atan(1.);
  double ynew, yold, multiplier ;
  int flag;

  /* Allocate NR vectors for maximum size used */
  T=dvector(0,N);
  y=dvector(1,2);
  dydx=dvector(1,2);
  
  //check if k1-k2 is legal

  if((abs(k1-k2)%2) != 0){
    printf("k1 - k2 is not an even number. Any self respecting 4th order Run\
ge-Kutta algorithm cannot proceed with an odd difference !!");
    exit(1);
  }


  if(k1>k2){
    multiplier=1.;
    for(i=k1;i>=k2;i--)
      T[i] = E - V[i];
  }else{
    multiplier = -1.;
    for(i=k1;i<=k2;i++)
        T[i]=E-V[i];
  }

  //initial values 
  y[1]= *Psi;
  y[2]= *Psip;
  flag=1;

  if(Psiout != NULL)
    {
      Psiout[k1] = *Psi;
      ynew = *Psi;
    }else{
      ynew = *Psi;
    }
  if(k1<k2){  
    int dk = 2;
    double h = 2./N;
    double x;
    
    for (k=k1; k<=k2 - 2; k+=dk) {
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480sch(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      if(Psiout!=NULL)
	{
	  Psiout[k+2] = y[1];
	  yold = ynew;
	  ynew = y[1];
	}else{
	  yold = ynew;
	  ynew = y[1];
	}
      if(fabs(y[1])>SQRTBIG){
	y[1]=y[1]*SMALL;
	y[2]=y[2]*SMALL;	
	yold=yold*SMALL;
        ynew=ynew*SMALL;

	if(Psiout != NULL)
	  {
	    for(g=k1;g<=k+2;g+=dk){
	      Psiout[g]=Psiout[g]*SMALL;
	    }
	  }
      }else{
	if(fabs(y[1])<SQRTSMALL){
	  y[1]=y[1]*BIG;
	  y[2]=y[2]*BIG;
	  yold=yold*BIG;
	  ynew=ynew*BIG;

	  if(Psiout != NULL){
	    for(g=k1;g<=k+2;g+=dk){
	      Psiout[g]=Psiout[g]*BIG;
	    }
	  }
	}
      }
      if (ynew==0){
        nodes++;
	flag = 0;
      }else{
        if(flag != 0 && (ynew/yold)<0){
          nodes++;
        }else{flag = 1;}
      }
    }
  }

  //  flag = 1;

  
  if(k1>k2){
    int dk = -2;
    double h = -2./N;
    double x;
    for (k=k1; k>=k2+2; k+=dk) {
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      if(Psiout != NULL){
	Psiout[k-2] = y[1];
	yold = ynew;
	ynew = y[1];
      }else{
	yold = ynew;
	ynew = y[1];
      }
      if(fabs(y[1])>SQRTBIG){
	y[1]=y[1]*SMALL;
	y[2]=y[2]*SMALL;
	yold=yold*SMALL;
	ynew=ynew*SMALL;
	if(Psiout != NULL){
	  for(g=k1;g>=k-2;g+=dk){
	    Psiout[g]=Psiout[g]*SMALL;
	  }
	}
      }else{
	if(fabs(y[1])<SQRTSMALL){
	  y[1]=y[1]*BIG;
	  y[2]=y[2]*BIG;
	  yold=yold*BIG;
	  ynew=ynew*BIG;
	  if(Psiout != NULL){
	    for(g=k1;g>=k-2;g+=dk){
	      Psiout[g]=Psiout[g]*BIG;
	    }
	  }
	}
      }
      if(ynew==0){
	nodes++;
	//flag = 0;
      }else{
	if(/*flag != 0 && */ (ynew/yold)<0){
	  nodes++;
	}//else{/*flag = 1;*/}
      }
    }
  }
  
  *Psi = y[1];
  *Psip = y[2];
  return nodes;
  
  free_dvector(T,0,N);
  free_dvector(y,1,2);
  free_dvector(dydx,1,2);

}

double func_Schrodinger(double E, int match, double V[], double r[], double dr[], int N)
{
  int k1, k2;
  double Psi, Psip;

  k1 = 0;
  k2 = N;
  Psi = 0.;
  Psip = 1.;
  schint(&Psi,&Psip,NULL,k1,k2,V,E,r,dr,N);
  return Psi;  
}

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N)
{
  int nodes;
  int k1, k2;
  double Psi, Psip;
  k1 = 0 ;
  k2 = N;
  Psi = 0.;
  Psip = 1.;
  nodes = schint(&Psi,&Psip,NULL,k1,k2,V,E,r,dr,N);
  return (double)(nodes-match);
}

#define TOL (1.e-12)

void getEs(double E[], int nmax, double Elower, double V[], double r[], double dr[], int N)
{
  int i;
  double E1,E2;
  for(i=0;i<=nmax;i++){
    E1=rtbisp480(func_SchrodingerNodes,Elower,0.,TOL,i,V,r,dr,N);
    E2=rtbisp480(func_SchrodingerNodes,Elower,0.,TOL,i+1,V,r,dr,N);
    E[i]=zriddrp480(func_Schrodinger,E1,E2,TOL,0,V,r,dr,N);
    Elower = E[i];
  }
}

void getallEs(double **E, int lmax, int nmax[], double Z, double V[], double r[], double dr[], int N)
{
  double *Veff;
  int k,lcount;
  double l;
  Veff = dvector(0,N);
  for(lcount=0;lcount<=lmax;lcount++){
    l = (double) lcount;
    for(k=1;k<=N;k++){
      Veff[k] = V[k] + l*(l+1)/(2*r[k]*r[k]);
    }
    Veff[0] = 0.;
    getEs(E[lcount],nmax[lcount],-Z*Z,Veff,r,dr,N);
  }
  free_dvector(Veff,0,N);
}

void getPsi(double Psi[], double E, int l, double V[], double r[], double dr[], int N)
{
  int k,count1,count2,count3,count4,count5;
  double *Veff ;
  double Psif,Psip,*Psimodsq;
  int kmatch=0;
  double normfactor;
  double ld;

  Veff = dvector(0,N);
  Psimodsq = dvector(0,N);

  ld = (double)l;
  for(k=0;k<=N;k++){
    Veff[k] = V[k] + ld*(ld+1)/(2*r[k]*r[k])  ;
  }
  Veff[0]=0.;

  for(count1=0;count1<=N;count1+=2)
    if(Veff[count1]<E) kmatch=count1;
  if(kmatch==0){
    printf("Veff never below E=%f in getPsi. \n",E);
    exit(1);
  }

  Psif = 0.;
  Psip = 1.;

  schint(&Psif,&Psip,Psi,0,kmatch,Veff,E,r,dr,N);
  for(count2=0;count2<=kmatch;count2+=2){
    Psi[count2] = Psi[count2]/Psi[kmatch];
  }

  Psif = 0.;
  Psip = 1.;
  schint(&Psif,&Psip,Psi,N,kmatch,Veff,E,r,dr,N);
  for(count3=N;count3>=kmatch;count3 -= 2){
    Psi[count3] = Psi[count3]/Psi[kmatch];
  }
  
  //interpolating the odd points
    interpolate(Psi,N);

  //calculating the probability at all r's
  for(count4=0;count4<=N;count4++){
  Psimodsq[count4] = Psi[count4]*Psi[count4];
  }
  
  //integrating the probability for normalisation
  normfactor = sqrt(simpint(Psimodsq,r,dr,N));

  //normalising the wavefunction
  for(count5=0;count5<=N;count5++){
    Psi[count5] = Psi[count5]/normfactor;
  }

  free_dvector(Veff,0,N);
  free_dvector(Psimodsq,0,N);
}

void getallPsis(double ***Psi, double **E,int lmax, int *nmax,double *V,double r[],double dr[], int N)
{
  int lc, n;
  double l;
  for(lc=0;lc<=lmax;lc++){
    l = (double)lc;
    for(n=0;n<=nmax[lc];n++)
      getPsi(Psi[lc][n],E[lc][n],l,V,r,dr,N);
  }
}

void getRho(double *Rho, double ***Psi, double **F, int lmax, int *nmax, int N){
  int k,l,n;
  
  for(k=0;k<=N;k++)
    {Rho[k]=0.;
      for(l=0;l<=lmax;l++)
	for(n=0;n<=nmax[l];n++)
	  Rho[k]=Rho[k] + F[l][n]*Psi[l][n][k]*Psi[l][n][k];
      
    }
}

void derivs_poisson(double x, double y[], double dydx[], int k, double Rho[], double r[], double dr[], int N)
{
  dydx[1] = y[2]*dr[k];
  dydx[2] = -Rho[k]*dr[k]/r[k];
}

void getphi(double phi[], double Rho[], double r[], double dr[], int N)
{
    
    /* Working variables */
    double *y,*dydx,*Phi; /* NR vector for diff eq's */
    double x,h;
    int k,dk;

    /* Value of pi */

    double pi;
    pi=4.*atan(1.);

    /* Allocate NR vectors for maximum size used */
    Phi=dvector(0,N);
    y=dvector(1,2);
    dydx=dvector(1,2);

      /* Runge-Kutta solution using rk4p480(). */
      /* Set up initial step sizes (h,dk) and initial conditions ... */
      h = -2./(double)N;
      dk = -2;
      y[1]=simpint(Rho,r,dr,N);
      y[2]=0.;
      Phi[N]=y[1];
      for (k=N; k>=2; k+=dk) {
	x=-k*h;
	derivs_poisson(x,y,dydx,k,Rho,r,dr,N);
	rk4p480sch(y,dydx,2,x,h,y,derivs_poisson,k,dk,Rho,r,dr,N);
	Phi[k+dk]=y[1];
      }
      
      interpolate(Phi,N);
      
      for(k=N; k>=0;k--){
	phi[k]=Phi[k]/r[k];
      }
      phi[0]=0.;
      
      /* Deallocate NR vectors: Should always clean up space at the end. */
      free_dvector(Phi,0,N);
      free_dvector(y,1,2);
      free_dvector(dydx,1,2);
}

double excPZ(double rs){
  double exc;
  double pi = 4.*atan(1.);
  if (rs<1.)
    {exc = -(3./4.)*pow((3./(2.*pi)),(2./3.))/rs + 0.0311*log(rs) + (-0.0480) + 0.002*rs*log(rs)+(-0.0116)*rs;
    }else{
      exc = -(3./4.)*pow((3./(2.*pi)),(2./3.))/rs + (-0.1423)/(1. + 1.0529*pow(rs,0.5) + 0.3334*rs);
    }
  return exc;
}

double excpPZ(double rs){
  double excp;
  double pi = 4.*atan(1.);
  if (rs<1.)
    {excp = (3./4.)*pow((3./(2.*pi)),(2./3.))/(rs*rs) + 0.0311/rs + 0.002*(1.+log(rs))+(-0.0116);
    }else{
      excp = (3./4.)*pow((3./(2.*pi)),(2./3.))/(rs*rs) - (-0.1423)/pow((1. + 1.0529*pow(rs,0.5) + 0.3334*rs),2.)*(1.0529/(2.*pow(rs,0.5))+0.3334);
    }
  return excp;
}

void getVxc(double Vxc[], double Rho[], double r[], double dr[], int N)
{ int k;
 double realrho, rs;
 double pi = 4.*atan(1.);
 for(k=1;k<=N;k++){
   if(fabs(Rho[k])<=SMALL){
     Vxc[k]=0.;}
   else{
     realrho=Rho[k]/(4.*pi*r[k]*r[k]);
     rs=pow((3./(4.*pi)),1./3.)*pow(realrho,-1./3.);
     Vxc[k]=excp(rs)*(-rs/3.) + exc(rs);
   }
 }
 Vxc[0]=0.;
}

void getDepsxc(double Depsxc[], double Rho[], double r[], double dr[], int N)
{ int k;
 double pi = 4.*atan(1.); 
 double realrho, rs;

 for(k=1;k<=N;k++){
   if(fabs(Rho[k])<=SMALL){
     Depsxc[k]=0.;
   }else{
     realrho=Rho[k]/(4.*pi*r[k]*r[k]); 
     rs=pow((3./(4.*pi)),1./3.)*pow(realrho,-1./3.);
     Depsxc[k]=excp(rs)*(rs/3.) ;
   }
 }
 Depsxc[0]=0.;
}

double exc(double rs)
{
  /* constants */
  const double
    pi = 4.*atan(1.),
    X1 = 0.75*pow(3.0/(2.0*pi),2.0/3.0),  /* Exchange energy coeff */
    A  =  0.0310907,
    x0 = -0.10498,
    b  = 3.72744,
    c  = 12.9352,
    Q  = sqrt(4*c-b*b),
    X0 = x0*x0+b*x0+c;

  double x=sqrt(rs),X=x*x+b*x+c;

  return -X1/rs 
    + A*(
	 +log(x*x/X)+2*b/Q*atan(Q/(2*x+b))
	 -(b*x0)/X0*(
		     log((x-x0)*(x-x0)/X)+2*(2*x0+b)/Q*atan(Q/(2*x+b))
		     )
	 );
}

double excp(double rs)
{
  /* constants */
  const double
    pi = 4.*atan(1.),
    X1 = 0.75*pow(3.0/(2.0*pi),2.0/3.0),  /* Exchange energy coeff */
    A  =  0.0310907,
    x0 = -0.10498,
    b  = 3.72744,
    c  = 12.9352,
    Q  = sqrt(4*c-b*b),
    X0 = x0*x0+b*x0+c;
  
  double x=sqrt(rs),X=x*x+b*x+c;
  
  double dx=0.5/x; /* Chain rule needs dx/drho! */
  
  return dx*(
	     2*X1/(rs*x)+A*(
			    2./x-(2*x+b)/X-4*b/(Q*Q+(2*x+b)*(2*x+b))
			    -(b*x0)/X0*(2/(x-x0)-(2*x+b)/X-4*(2*x0+b)/
					(Q*Q+(2*x+b)*(2*x+b)) )
			    )
	     );
}

double getg(double g[], double Rho[],double Z, int lmax, int nmax[], int nmaxmax, double **F, double r[], double dr[], int N)
{
  int k,l,n;
  double *phi, *Vxc, ***Psi, *Rhonew, *V, **E, *Depsxc, *integrand;
  double Etot;

  phi=dvector(0,N);
  V=dvector(0,N);
  Vxc=dvector(0,N);
  Rhonew=dvector(0,N);
  Psi=d3tensor(0,lmax,0,nmaxmax,0,N);
  E=dmatrix(0,lmax,0,nmaxmax);
  Depsxc=dvector(0,N);
  integrand=dvector(0,N);

  getphi(phi,Rho,r,dr,N);
  getVxc(Vxc,Rho,r,dr,N);
  for (k=0; k<=N; k++)
    V[k]=-Z/r[k]+phi[k]+Vxc[k];
  V[0]=0.;

  getallEs(E,lmax,nmax,Z,V,r,dr,N);
  getallPsis(Psi,E,lmax,nmax,V,r,dr,N);
  getRho(Rhonew,Psi,F,lmax,nmax,N);

  for(k=0;k<=N;k++) g[k]=Rhonew[k]-Rho[k];

  getDepsxc(Depsxc,Rho,r,dr,N);
  for (k=0; k<=N; k++)
    integrand[k]=(-0.5*phi[k]+Depsxc[k])*Rho[k];
  Etot=simpint(integrand,r,dr,N);
  
  /* Add on the sum of the electron energies times the occupancies */
  for (l=0; l<=lmax; l++)
    for (n=0; n<=nmax[l]; n++)
      Etot+=F[l][n]*E[l][n];
  
  return Etot;

  free_dvector(V,0,N);
  free_dmatrix(E,0,lmax,0,nmaxmax);
  free_d3tensor(Psi,0,lmax,0,nmaxmax,0,N);
  free_dvector(Rhonew,0,N);
  free_dvector(phi,0,N);
  free_dvector(Vxc,0,N);
  free_dvector(Depsxc,0,N);
  free_dvector(integrand,0,N);
 
}
