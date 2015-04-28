#include <stdio.h>
#include <math.h>
#include "p480.h"
#include "nrutil.h"

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

/* ============ end of function declarations================================ */



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
  double ynew, yold ;
  int flag;

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
      printf("k1 - k2 is not an even number. Any self respecting 4th order Runge-Kutta algorithm cannot proceed with an odd difference !!");
      exit(1);
    }
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
      rk4p480(y,dydx,2,x,h,y,derivs_Schrodinger,k,dk,T,r,dr,N);
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
	printf("announcing the explosion of y!!\n");
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
	  printf("announcing the collapse of y!!");
	  y[1]=y[1]*BIG;
	  y[2]=y[2]*BIG;
	  yold=yold*BIG;
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
  //double *Psiout;
  int nodes;
  //Psiout = dvector(0,N);
  int k1, k2;
  double *Psi, *Psip;
  //Psip = (double *)malloc(sizeof(double));
  Psi = dvector(0,0);
  Psip = dvector(0,0);
  k1 = 0.;
  k2 = N;
  *Psi = 0.;
  *Psip = 1.;
  nodes = schint(Psi,Psip,NULL,k1,k2,V,E,r,dr,N);
  return *Psi;  
  //free_dvector(Psiout,0,N);
  free_dvector(Psi,0,0);
  free_dvector(Psip,0,0);
}

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N)
{
  //double *Psiout;
  int nodes;
  //Psiout = dvector(0,N);
  int k1, k2;
  double *Psi, *Psip;
  Psi = dvector(0,0);
  Psip = dvector(0,0);
  k1 = 0.;
  k2 = N;
  *Psi = 0.;
  *Psip = 1.;
  nodes = schint(Psi,Psip,NULL,k1,k2,V,E,r,dr,N);
  return (double)(nodes-match);
  //  free_dvector(Psiout,0,N);
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
    //    printf("done with rtb and zriddr\n");
    Elower = E[i];
    //printf("this is the last st of getEs\n");
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
  double *Psif,*Psip,*Psimodsq;
  int kmatch=0;
  double normfactor;
  double ld;

  Veff = dvector(0,N);
  Psimodsq = dvector(0,N);
  Psif = dvector(0,0);
  Psip = dvector(0,0);

  ld = (double)l;
  printf("l=%d\n",l);
  for(k=0;k<=N;k++){
    Veff[k] = V[k] + ld*(ld+1)/(2*r[k]*r[k])  ;
  }
  Veff[0]=0.;
  //  for(k=0;k<=N;k++) printf("%f\t %f\n",V[k],Veff[k]);
  for(count1=0;count1<=N;count1+=2)
    if(Veff[count1]<E) kmatch=count1;
  if(kmatch==0){
    printf("Veff never below E=%f in getPsi. \n",E);
    exit(1);
  }

  printf("made kmatch = %d\n",kmatch);
  *Psif = 0.;
  *Psip = 1.;

  schint(Psif,Psip,Psi,0,kmatch,Veff,E,r,dr,N);
  for(count2=0;count2<=kmatch;count2+=2){
    Psi[count2] = Psi[count2]/Psi[kmatch];
  }

  *Psif = 0.;
  *Psip = 1.;
  schint(Psif,Psip,Psi,N,kmatch,Veff,E,r,dr,N);
  for(count3=N;count3>=kmatch;count3 -= 2){
    Psi[count3] = Psi[count3]/Psi[kmatch];
  }
  
  //interpolating the odd points
  interpolate(Psi,N);
  //  for(k=0;k<=N;k++)
    //  printf("%f\n",Psi[k]);
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

  free_dvector(Psif,0,0);
  free_dvector(Psip,0,0);
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
      for(l=0;l<=lmax;l++){
	for(n=0;n<=nmax[l];n++){
	  Rho[k]=Rho[k] + F[n][l]*Psi[l][n][k]*Psi[l][n][k];
	}
      }
    }
}

void derivs_poisson(double x, double y[], double dydx[], int k, double Rho[], double r[], double dr[], int N)
{
  double pi;
  pi = 4.*atan(1.);
  dydx[1] = y[2]*dr[k];
  dydx[2] = -Rho[k]*dr[k]/r[k];
}

void getphi(double phi[], double Rho[], double r[], double dr[], int N)
{
    
    /* Working variables */
    int i;
    double *y,*dydx,*Phi,*integrand; /* NR vector for diff eq's */
    double x,h;
    int k,dk;

    /* Value of pi */

    double pi;
    pi=4.*atan(1.);

    /* Allocate NR vectors for maximum size used */
    Phi=dvector(0,N);
    y=dvector(1,2);
    dydx=dvector(1,2);
    integrand=dvector(0,N);

      /* Runge-Kutta solution using rk4p480(). */
      /* Set up initial step sizes (h,dk) and initial conditions ... */
      h = -2./(double)N;
      dk = -2;
      for(k=0;k<=N;k++)integrand[k]=4*pi*Rho[k]*r[k]*r[k];
      y[1]=1.;
      y[2]=0.;
      Phi[N]=1;
      for (k=N; k>=2; k+=dk) {
	x=-k*h;
	derivs_poisson(x,y,dydx,k,Rho,r,dr,N);
	rk4p480(y,dydx,2,x,h,y,derivs_poisson,k,dk,Rho,r,dr,N);
	Phi[k-2]=y[1];
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
    free_dvector(integrand,0,N);
}
