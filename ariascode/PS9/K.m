% Overlap operator (acting on 3d data sets)
%
% Usage: out=K(W)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_R: Lattice vectors

function out=K(W)
  
  %# Must declare all globals with such statements to access them
  global gbl_G2;

out=((ones(1,size(gbl_G2,1))'./(gbl_G2+ones(1,size(gbl_G2,1))'))*ones(1,size(W,2))).*W;%# <=== YOUR CODE HERE

endfunction


