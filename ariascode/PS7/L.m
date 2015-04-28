% Overlap operator (acting on 3d data sets)
%
% Usage: out=L(in)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_R: Lattice vectors

function out=L(in)
  global gbl_R; %# Must declare all globals with such statements to access them
  global gbl_G2;

  %# Operator definition (multiplication by volume)
out=-det(gbl_R)*gbl_G2.*in;%# <=== YOUR CODE HERE
endfunction


