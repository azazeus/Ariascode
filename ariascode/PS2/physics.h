int schint(double *Psi, double *Psip, double Psiout[], int k1, int k2, double V[], double E, double r[], double dr[], int N);

double func_Schrodinger(double E, int match, double V[], double r[], double dr[], int N);

double func_SchrodingerNodes(double E, int match, double V[], double r[], double dr[], int N);

void getES(double E[], int nmax, double Elower, double V[], double r[], double dr[], int N);
