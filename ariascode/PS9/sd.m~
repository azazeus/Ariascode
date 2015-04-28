%# method of steepest descent for finding minima of Energy

function out = sd(W,Nit)

alpha = 3e-5;
for i=1:Nit;
  W=W-alpha*getgrad(W);
  out=W;
%#printf(' Nit = %d \t E = %20.16f \n',i, getE(W));
end

endfunction