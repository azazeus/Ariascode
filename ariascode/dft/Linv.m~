% Overlap operator (acting on 3d data sets)
%
% Usage: out=Linv(in)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_R: Lattice vectors
% gbl_G2: sq. magnitude of reciprocal lattice vectors
% gbl_S
function out=Linv(in)
  global gbl_R; %# Must declare all globals with such statements to access them
  global gbl_G2;
  global gbl_S;

  %# Operator definition (multiplication by volume)
out=-1/det(gbl_R)*(ones(prod(gbl_S),1)./gbl_G2).*in;%# <=== YOUR CODE HERE
out(1)=0;
endfunction


