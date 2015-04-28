% FFT operator (acting on 3d data sets)
%
% Usage: out=cI(in)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_S
function out=cIdag(in)
  global gbl_S;
  out = zeros(size(in));

  %# Operator definition (multiplication by volume)

for col=1:size(in,2)
  out(:,col)=fftw3(in(:,col),gbl_S(1),gbl_S(2),gbl_S(3),-1);%# Columnwise
end

endfunction


