% IFFT operator (acting on 3d data sets)
%
% Usage: out=cJ(in)
%
% in: input 3d data set
% out: output 3d data set
%
% Uses GLOBAL variable(s) ---
% gbl_S
function out=cJdag(in)
  global gbl_S;
  out=zeros(size(in));
  %# Operator definition (multiplication by volume)

for col=1:size(in,2)
out(:,col)=1/(gbl_S(1)*gbl_S(2)*gbl_S(3))*fftw3(in(:,col),gbl_S(1),gbl_S(2),gbl_S(3),1);
end

endfunction


