%# method of steepest descent for finding minima of Energy using getE()
 # and getgrad()
%# usage : [out, Elist] = sd(W, Nit)
function [out, Elist] = sd(W,Nit)

alpha = 3e-5;
out=W;
for i=1:Nit;
  out=out-alpha*getgrad(out);
Elist(i)=getE(out);
%#printf(' Nit = %d \t E = %20.16f \n',i, getE(W));
Elist(i)
end

endfunction