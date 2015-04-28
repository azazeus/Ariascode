%# program to perform line minimisation 
%# 
%# usage [out, Elist]=lm(W, Nit)

function [out, Elist]=lm(W, Nit)
out=W;
alphat=9e-5;
d=W; %# filler for declaring d
for it=1:Nit
g=getgrad(out);
if(it>1)
 printf('angle cosine : %f and energy : %f \
	\n',real(sum(diag(g'*d)))/(sqrt(real(sum(diag(d'*d)))*real(sum(diag(g'*g))))), \
	getE(out));
endif
d=-g;
gt=getgrad(out+alphat*d);
alpha=alphat*(real(sum(diag(g'*d))))/(real(sum(diag((g-gt)'*d))));
out=out+alpha*d;
Elist(it)=getE(out);
end

endfunction
