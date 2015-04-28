#ifndef _PHYSICS_H_
#define _PHYSICS_H_
void derivs_poisson(double x,double y[],double dydx[],int k,double rho[],double r[],double dr[],int N);

void getphi(double phi[],double Rho[],double r[],double dr[],int N);

void getRho(double *Rho,double ***Psi,double **F,int lmax,int *nmax,int N);

void getallPsis(double ***Psi,double **E,int lmax,int *nmax,double *V,double r[],double dr[],int N);

void getPsi(double Psi[], double E, int l, double V[], double r[], double dr[], int N);

void getallEs(double **E, int lmax, int nmax[], double Z, double V[], double r[], double dr[], int N);

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N);

void getEs(double E[], int nmax, double Elower, double V[], double r[], double dr[], int N);

double func_Schrodinger(double E, int match, double V[], double r[], double dr[], int N);

void derivs_Schrodinger(double x, double y[], double dydx[], int k, double T[], double r[], double dr[], int N);

int schint(double *Psi, double *Psip, double Psiout[],int k1, int k2, double V[], double E,double r[], double dr[], int N);

double excPZ(double rs);

double excpPZ(double rs);

void getVxc(double Vxc[], double Rho[], double r[], double dr[], int N);

void getDepsxc(double Depsxc[], double Rho[], double r[], double dr[], int N);

double exc(double rs);

double excp(double rs);

double getg(double g[], double Rho[], double Z, int lmax, int nmax[], int nmaxmax, double **F, double r[], double dr[], int N);

#endif
