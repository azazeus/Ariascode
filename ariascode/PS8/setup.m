%# Make setup info globally accessible (ugh!!#$@#%)
global gbl_S; global gbl_R; global gbl_G2;
S=[20; 25; 30]
R=diag([6;6;6]);
%# Code fragment to create columns of m1, m2, m3 indices and the matrix
 # M
ms=[0:prod(S)-1]';
m1=rem(ms,S(1));
m2=rem(floor(ms/S(1)),S(2));
m3=rem(floor(ms/(S(1)*S(2))),S(3));
M=[m1, m2, m3];
n1=m1-(m1>S(1)/2)*S(1);
n2=m2-(m2>S(2)/2)*S(2);
n3=m3-(m3>S(3)/2)*S(3);
N=[n1,n2,n3];
r=M*inv(diag(S))*R';
G=2*pi*N*inv(R)';
G1=G.^2;
G2=sum(G1,2);
%# Assign computed values to the global variables
gbl_S=S; gbl_R=R; gbl_G2=G2;