%# this function computes the potential V for harmonic oscillator with 
				# omega = 2

global gbl_Vdual;
global gbl_f;
global gbl_R;
global gbl_G;
gbl_f=2 %# the usual case of constant filling of 2 electrons per orbital

y=1.75;

ZI=1.0;
rI=[8;8;8];
rI2=[8;8;8]+[y;0;0];

gbl_Vdual=cJdag(O(Linv(4*pi*O(ZI/det(gbl_R)*(exp(-i*gbl_G*rI)+exp(-i*gbl_G*rI2))))));

%# Finite difference test

Ns = 1 %# Number of states

randn('seed', 0.2004);
W=(randn(prod(S),Ns)+i*randn(prod(S),Ns));

more off; %# View output as it is computed
fdtest(W);

%# orthogonalising W
N=W'*O(W);
TempW=W*sqrtm(inv(N));
W=TempW;

%# Converge using steepest descent
%# W=lm(W,20);
%# [Wlm, Elm]=lm(W,50);
%# [Wpclm, Epclm]=pclm(W,50);
%# [Wcg1, Ecg1]=pccg(W,50,1);
%#[Wcg2, Ecg2]=pccg(W,50,2);
%#[Wcg3, Ecg3]=pccg(W,50,3);
%#its=[1:50];
%#semilogy(its,Elm-18,its,Epclm-18,its,Ecg1-18,its,Ecg2-18,its,Ecg3-18);
%# Set W to best result for plotting
%#W=Wcg3;
%#pause;

%# lots of digits
format long

%# converge a little using steepest descents
W=sd(W,10);

%# converge more with preconditioned line minimisation
W=pclm(W,10);

%# finish convergence with full preconditioned conjugate gradients
W=pccg(W,15,3);

%# Extract and display final results
[Psi, epsilon] = getPsi(W);

for st=1:Ns
  printf('===State # %d, Energy = %f ==== \n', st, epsilon(st));
	%# viewmid(abs(cI(Psi(:,st))).^2,S)
end