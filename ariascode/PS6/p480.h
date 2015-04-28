#ifndef _P480_H_
#define _P480_H_

double trapint(double f[], double r[], double dr[], int N);
double simpint(double f[], double r[], double dr[], int N);

void rk4p480(double y[], double dydx[], int n, double x, double h, double yout[],
	     void (*derivs)(double, double [], double [],
			    int, double [], double [], double [], int),
	     int k, int dk, double f[], double r[], double dr[], int N);

double rtbisp480(double (*func)(double,int,double [],double [], double [], int),
		 double x1, double x2, double xacc, 
		 int match, double V[], double r[], double dr[], int N);

double zriddrp480(double (*func)(double,int,double [],double [],double [], int),
		  double x1, double x2, double xacc,
		  int match, double V[], double r[], double dr[], int N);

void interpolate(double f[], int N);

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);

void ludcmpp480(double **a, int n, int *indx, double *d);

void lubksbp480(double **a, int n, int *indx, double b[]);


#endif

