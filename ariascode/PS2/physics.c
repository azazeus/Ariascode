#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

int schint(double *,double *,double *,int,int,double *,double,double *,double *,int);

double func_Schrodinger(double, int, double *, double *, double *, int);

double func_SchrodingerNodes(double, int, double *, double *, double *, int);

void getES(double *, int, double, double *, double *, double *, int);

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
    
  /* Allocate NR vectors for maximum size used */
  T=dvector(0,N);
  y=dvector(1,2);
  dydx=dvector(1,2);
  
  //check if k1-k2 is legal
  if(k1>k2){
    for(i=k1;i>=k2;i--)
      {  T[i] = E - V[i];
      }
    if(((k1-k2)%2) != 0){
      printf("k1 - k2 is not an even number. Any self respecting 4th order Runge-Kutta algorithm cannot proceed with an odd difference !!");
      exit(1);
    }
  }else{
    for(i=k1;i<=k2;i++)
      {  T[i]=E-V[i];
      }
    if(((k2-k1)%2)!=0){
      printf("k1 - k2 is not an even number. Any self respecting 4th order Rung\e-Kutta algorithm cannot proceed with an odd difference !!");
      exit(1);
    }
  }

  //initial values 
  y[1]= Psiout[k1] = *Psi;
  y[2]= *Psip;
  
  if(k1<k2){  
    int dk = 2;
    double h = 2./N;
    double x;
    
    for (k=k1; k<=k2 - 2; k+=dk) {
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      Psiout[k+2] = y[1];
      if(fabs(y[1])>SQRTBIG){
	y[1]=y[1]*SMALL;
	y[2]=y[2]*SMALL;	
	for(g=k1;g<=k+2;g+=dk){
	  Psiout[g]=Psiout[g]*SMALL;
	}
      }else{
	if(fabs(y[1])<SQRTSMALL){
	y[1]=y[1]*BIG;
	y[2]=y[2]*BIG;
	  for(g=k1;g<=k+2;g+=dk){
	    Psiout[g]=Psiout[g]*BIG;
	  }
	}
      }
      if (Psiout[k+2]==0){
        nodes++;
      }else{
        if((Psiout[k+2]/Psiout[k])<0){
          nodes++;
        }
      }
    }
  }
  
  if(k1>k2){
    int dk = -2;
    double h = -2./N;
    double x;
    for (k=k1; k>=k2+2; k+=dk) {
      x=k*h;
      derivs_Schrodinger(x,y,dydx,k,T,r,dr,N);
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
      Psiout[k-2] = y[1];
      if(fabs(y[1])>SQRTBIG){
	printf("announcing the explosion of y!!");
	y[1]=y[1]*SMALL;
	y[2]=y[2]*SMALL;
	for(g=k1;g<=k-2;g+=dk){
	  Psiout[g]=Psiout[g]*SMALL;
	}
      }else{
	if(fabs(y[1])<SQRTSMALL){
	  printf("announcing the collapse of y!!");
	  y[1]=y[1]*BIG;
	  y[2]=y[2]*BIG;
	  for(g=k1;g<=k-2;g+=dk){
	    Psiout[g]=Psiout[g]*BIG;
	  }
	}
      }
      if (Psiout[k-2]==0){
	nodes++;
      }else{
	if((Psiout[k-2]/Psiout[k])<0){
	  nodes++;
	}
      }
    }
  }

  *Psi = y[1];
  *Psip = y[2];
  return nodes;
  
  free_dvector(T,0,N);
  free_dvector(y,1,1);
  free_dvector(dydx,1,1);

}

double func_Schrodinger(double E, int match, double V[], double r[], double dr[], int N)
{
  double *Psiout;
  int nodes;
  Psiout = dvector(0,N);
  int k1, k2;
  double *Psi, *Psip;
  //Psip = (double *)malloc(sizeof(double));
  Psi = dvector(0,0);
  Psip = dvector(0,0);
  k1 = 0.;
  k2 = N;
  *Psi = 0.;
  *Psip = 1.;
  nodes = schint(Psi,Psip,Psiout,k1,k2,V,E,r,dr,N);
    return *Psi;  
  free_dvector(Psiout,0,N);
  free_dvector(Psi,0,0);
  free_dvector(Psip,0,0);
}

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N)
{
  double *Psiout;
  int nodes;
  Psiout = dvector(0,N);
  int k1, k2;
  double *Psi, *Psip;
  Psi = dvector(0,0);
  Psip = dvector(0,0);
  k1 = 0.;
  k2 = N;
  *Psi = 0.;
  *Psip = 1.;
  nodes = schint(Psi,Psip,Psiout,k1,k2,V,E,r,dr,N);
  return (double)(nodes-match);
  free_dvector(Psiout,0,N);
  free_dvector(Psi,0,0);
  free_dvector(Psip,0,0);
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
