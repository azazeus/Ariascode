%# this function computes the potential V for harmonic oscillator with 
				# omega = 2

global gbl_Vdual;


V=2*gbl_dr.^2;
gbl_Vdual=cJdag(O(cJ(V)));

%# Finite difference test

Ns = 4 %# Number of states

randn('seed', 0.2004);
W=(randn(prod(S),Ns)+i*randn(prod(S),Ns));

more off; %# View output as it is computed
fdtest(W);

%# orthogonalising W
N=W'*O(W);
TempW=W*sqrtm(inv(N));
W=TempW;

%# Converge using steepest descent
W=sd(W,400);

%# Extract and display final results
[Psi, epsilon] = getPsi(W);

for st=1:Ns
  printf('===State # %d, Energy = %f ==== \n', st, epsilon(st));
	 viewmid(abs(cI(Psi(:,st))).^2,S)
end