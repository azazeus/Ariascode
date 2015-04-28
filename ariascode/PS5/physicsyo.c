#include <stdio.h>
#include "nrutil.h"
#include <math.h>
#include "p480.h"
#include "physics.h"

#define TOL (1.e-12)
#define SMALL (1.e-300)
#define SQRTSMALL (1.e-150)
#define BIG (1.e300)
#define SQRTBIG (1.e150)

void derivs_poisson(double x,double y[],double dydx[],int k,double Rho[],double r[],double dr[],int N) {
  dydx[1]=y[2]*dr[k];
  dydx[2]=-Rho[k]/r[k]*dr[k];
}

void getphi(double phi[],double Rho[],double r[],double dr[],int N){
  int k,dk;
  double x,h;
  double *y,*dydx;

  y=dvector(1,2);
  dydx=dvector(1,2);

  h=-2./N;
  dk=-2;
  y[1]=simpint(Rho,r,dr,N);
  y[2]=0.;
  phi[N]=y[1];
  for (k=N;k>=2;k+=dk){
    x=-k*h;
    derivs_poisson(x,y,dydx,k,Rho,r,dr,N);
    rk4p480(y,dydx,2,x,h,y,derivs_poisson,k,dk,Rho,r,dr,N);
    phi[k+dk]=y[1];
  }
    interpolate(phi,N);
    for (k=0;k<=N;k++)
      phi[k]/=r[k];
    phi[0]=0.;
    free_dvector(y,1,2);
    free_dvector(dydx,1,2);
}
    
void getRho(double *Rho,double ***Psi,double **F,int lmax,int *nmax,int N){
  int l,n,k;
  for (k=0;k<=N;k++)
    Rho[k]=0.;
  for (k=0;k<=N;k++)
    for (l=0;l<=lmax;l++)
      for (n=0;n<=nmax[l];n++)
	Rho[k]=Rho[k]+F[l][n]*Psi[l][n][k]*Psi[l][n][k];
}

void getallPsis(double ***Psi,double **E,int lmax,int *nmax, double *V,double r[],double dr[],int N){
  int l,n;
  for (l=0;l<=lmax;l++)
    for (n=0;n<=nmax[l];n++)
      getPsi(Psi[l][n],E[l][n],l,V,r,dr,N);
}

void getPsi(double Psi[], double E, int l, double V[], double r[], double dr[], int N) {
  double *Veff;
  int k;

  Veff=dvector(0,N);

  for (k=0;k<=N;k++) {
    Veff[k]=V[k]+l*(l+1)/(2*r[k]*r[k]);
  }
  Veff[0]=0.;

  int kmatch=0;
  for (k=0;k<=N;k++)
    if (Veff[k]<E) kmatch=k; 
  if (kmatch==0) {
    printf("Veff never below E=%f in getPsi. \n",E);
    exit(1);
  }
  if (kmatch%2 !=0) 
    kmatch--;
  double Psi_init=0.;
  double Psip_init=1.;
  double temp;
  schint(&Psi_init, &Psip_init, Psi, 0,kmatch, Veff,E,r,dr,N);
  temp=1./Psi[kmatch];
  for (k=0;k<=kmatch;k+=2)
    Psi[k]*=temp;
  Psi_init=0.;
  Psip_init=1.;
  schint(&Psi_init,&Psip_init, Psi,N,kmatch,Veff,E,r,dr,N);
  temp=1./Psi[kmatch];
  for (k=kmatch;k<=N;k+=2)
    Psi[k]*=temp;
  interpolate(Psi,N);

  double *Psi_squared;
  Psi_squared=dvector(0,N);
  for (k=0;k<=N;k++)
    Psi_squared[k]=Psi[k]*Psi[k];

  temp=simpint(Psi_squared,r,dr,N);
  temp=1/sqrt(temp);
  for (k=0;k<=N;k++)
    Psi[k]*=temp;
  
  free_dvector(Veff,0,N);
  free_dvector(Psi_squared,0,N);
}

void getallEs(double **E, int lmax, int nmax[], double Z, double V[], double r[], double dr[], int N) {
  int k,l;
  double *Veff;
  int nmaxmax=0;
  for (k=0;k<=lmax;k++)
    if (nmax[k] > nmaxmax) nmaxmax = nmax[k];
  Veff=dvector(0,N);
  for (l=0; l<=lmax;l++){
  for (k=0; k<=N;k++) {
    Veff[k]=V[k]+l*(l+1)/(2*r[k]*r[k]);
  }
  Veff[0]=0.;
  
  getEs(E[l],nmax[l],-Z*Z,Veff,r,dr,N);
  }
  free_dvector(Veff,0,N);
}

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N) {
  int k1;
  int k2;
  double Psi;
  double Psip;
  double *Psiout;
  double zero;

  Psiout=dvector(0,N);
  Psi=0.;
  Psip=1.;
  

  zero = (double)(schint(&Psi,&Psip,Psiout,k1=0,k2=N,V,E,r,dr,N));
  free_dvector(Psiout,0,N);
  zero =zero-match;
  return zero;

}

void getEs(double E[], int nmax, double Elower, double V[], double r[], double dr[], int N){
  int n;
  double E1=0.;
  double E2=0.;
  
  /* loop to get states */
  for (n=0; n<=nmax;n++) {
    /* Get E1 as an energy with n nodes */
    E1=rtbisp480(func_SchrodingerNodes,Elower,0.,TOL,n,V,r,dr,N);
    /* Get E2 as an energy with n+1 nodes */
    E2=rtbisp480(func_SchrodingerNodes,Elower,0.,TOL,n+1,V,r,dr,N);
    
    /* Now, get the solution, which is in between! */
    E[n]=zriddrp480(func_Schrodinger,E1,E2,TOL,0,V,r,dr,N);
  }
  
}

  

  
double func_Schrodinger(double E, int match, double V[], double r[], double dr[], int N) {
  int k1;
  int k2;
  double Psi;
  double Psip;
  double *Psiout;
  
  Psiout=dvector(0,N);
  Psi=0.;
  Psip=1.;
  

  schint(&Psi,&Psip,NULL,k1=0,k2=N,V,E,r,dr,N);
  free_dvector(Psiout,0,N);
  return Psi;

}


void derivs_Schrodinger(double x, double y[], double dydx[], int k, double T[], double r[], double dr[], int N) {
  dydx[1]=y[2] * dr[k];
  dydx[2]= -2. * T[k] * y[1] * dr[k];
}

int schint(double *Psi, double *Psip, double Psiout[], int k1, int k2, double V[], double E, double r[], double dr[], int N) {

  /* working variables */
  double *T;
  double *y,*dydx; /* NR vector for diff eq's */
  double x;
  double h;
  int k,dk,kk;
  int nodes = 0;

  if ((k2-k1)%2 != 0) {
    printf("Error in Schint: k2-k1 must be even in routine but k1=%d; and k2=%d\n",k1,k2);
    exit(1);
  }

  /* Allocate NR vectors for maximum size used */
  T=dvector(0,N);
  y=dvector(1,2);
  dydx=dvector(1,2);
  
  y[1] = *Psi;
  y[2] = *Psip;

  double Psiold;
  Psiold=y[1];
  if (k2 > k1) {
    h=2./N;
    dk=2;
    if (Psiout != NULL) Psiout[k1]=y[1];
    for (k=k1; k<=k2; k++){
      T[k]=E-V[k];
    }
    for (k=k1; k<=k2-2; k+=dk) {
      if (fabs(y[1]) > SQRTBIG) {
	y[1] = SMALL * y[1];
	y[2] = SMALL * y[2];
	if (Psiout != NULL)
	for (kk=k1;kk <= k; kk +=dk)
	  Psiout[kk]=SMALL*Psiout[kk];
      }
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      if (Psiout != NULL) Psiout[k+dk]=y[1];
      if (y[1] * Psiold < 0.) nodes++;
      Psiold=y[1];
    }
  }
  else {
    h=-2./N;
    dk=-2;
    if (Psiout != NULL) Psiout[k1]=y[1];
    for(k=k1; k>=k2; k--) {
      T[k]=E-V[k];
    }
    for (k=k1; k>=k2+2; k+=dk) {
      if (fabs(y[1]) > SQRTBIG){
	y[1]=SMALL * y[1];
	y[2]=SMALL * y[2];
	if (Psiout != NULL)
	for (kk=k1;kk >=k;kk +=dk)
	  Psiout[kk]=SMALL*Psiout[kk];
      }
      x=-h*k;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      if (Psiout != NULL) Psiout[k+dk]=y[1];
      if (y[1] * Psiold < 0.) nodes++;
      Psiold = y[1];
    }
  }
  *Psi = y[1];
  *Psip = y[2];
    
  /* Deallocate Nr vectors: Should always clean up space at the end. */
  free_dvector(T,0,N);
  free_dvector(y,1,2);
  free_dvector(dydx,1,2);
  return nodes;
}


double excPZ(double rs) {

  double pi=4.*atan(1.);
  double alpha=0.75*pow(3.0/(2.0*pi),2.0/3.0);
  double a = 0.0311;
  double b = -0.0480;
  double c = 0.0020;
  double d = -0.0116;
  double A = -0.1423;
  double B = 1.0529;
  double C = 0.3334;

  if ( rs < 1 )
    return (-alpha/rs + a*log(rs)+b+c*rs*log(rs)+d*rs);
  else
    return (-alpha/rs + A / (1. + B*sqrt(rs) + C*rs));
}

double excpPZ(double rs) {

  double pi=4.*atan(1.);
  double alpha=0.75*pow(3.0/(2.0*pi),2.0/3.0);
  double a = 0.0311;
  double b = -0.0480;
  double c = 0.0020;
  double d = -0.0116;
  double A = -0.1423;
  double B = 1.0529;
  double C = 0.3334;

  if (rs < 1)
    return (c + d + alpha/(rs*rs) + a/rs + c*log(rs));
  else
    return (alpha/(rs*rs) - A*(C + B/(2*sqrt(rs)))/((1+B*sqrt(rs)+C*rs)*(1+B*sqrt(rs)+C*rs)));

}

void getVxc(double Vxc[], double Rho[], double r[], double dr[], int N) {

  double pi=4.*atan(1.);
  int k;
  double n,rs,exc_value,excp_value;

  Vxc[0]=0.;

  for (k=1; k<=N; k++) {
    if (fabs(Rho[k]) > SMALL) {
      n = Rho[k]/(4.* pi * r[k] * r[k]);
      rs = pow(3/(4*pi),1./3.) * pow(n , -1./3.);
      exc_value = exc(rs);
      excp_value = excp(rs);
      Vxc[k] = exc_value - (1./3.) * excp_value * rs;
    }
    else Vxc[k] = 0.;
  }


}

void getDepsxc(double Depsxc[], double Rho[], double r[], double dr[], int N) {
  double pi=4.*atan(1.);
  int k;
  double n,rs,exc_value,excp_value;

  Depsxc[0]=0.;

  for (k=1; k<=N; k++) {
    if (fabs(Rho[k]) > SMALL) {
      n = Rho[k]/(4.* pi * r[k] * r[k]);
      rs = pow(3/(4*pi),1./3) * pow(n , -1./3.);
      excp_value = excp(rs);
      Depsxc[k] = (1./3.) * excp_value * rs;
    }
    else Depsxc[k] = 0.;
  }

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

double getg(double g[], double Rho[], double Z, int lmax, int nmax[], int nmaxmax, double **F, double r[], double dr[], int N) {

    

  /* Physics variables */
  double *V, *Rhonew, *Phi, *Vxc,*Depsxc, *integrand;
  double **E,***Psi;
  double Etot;

  /* Working variables */
  int n,l,k;
  double x;

  /* Solver iteration varibles */
  int it;

  /* Value of pi */
  const double pi=4.*atan(1.);


  /* The rest is now general for ANY case */
  E=dmatrix(0,lmax,0,nmaxmax); /* Make space for E's and Psi's */
  Psi=d3tensor(0,lmax,0,nmaxmax,0,N);
  Rhonew=dvector(0,N);
  Phi=dvector(0,N);
  Vxc=dvector(0,N);
  Depsxc=dvector(0,N);

  integrand=dvector(0,N);

  V=dvector(0,N);

      /* Make potential from zero charge (debugs NaNs, etc.) */
      getphi(Phi,Rho,r,dr,N);
      getVxc(Vxc,Rho,r,dr,N);
      for (k=0; k<=N; k++) {
        V[k]=-Z/r[k]+Phi[k]+Vxc[k];
      }
      V[0]=0.;

      getallEs(E,lmax,nmax,Z,V,r,dr,N);
      getallPsis(Psi,E,lmax,nmax,V,r,dr,N);
      getRho(Rhonew,Psi,F,lmax,nmax,N);

      /* Compute and output total energy */

      /* Get Correction to sum of electron energies */
      getDepsxc(Depsxc,Rho,r,dr,N);
      for (k=0; k<=N; k++)
        integrand[k]=(-0.5*Phi[k]+Depsxc[k])*Rho[k];
      Etot=simpint(integrand,r,dr,N);

      /* Add on the sum of the electron energies times the occupancies */
      for (l=0; l<=lmax; l++)
        for (n=0; n<=nmax[l]; n++)
          Etot+=F[l][n]*E[l][n];

      for (k=0;k<=N;k++)
      g[k]=Rhonew[k]-Rho[k];

    return Etot;
}
