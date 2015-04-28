% FFT operator (acting on 3d data sets)
%
% Usage: out=cI(in)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_S
function out=cI(in)
  global gbl_S;

  %# Operator definition (multiplication by volume)
out=fftw3(in,gbl_S(1),gbl_S(2),gbl_S(3),1);%# <=== YOUR CODE HERE

endfunction


