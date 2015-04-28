%# program to perform line minimisation 
%# 
%# usage [out, Elist]=pccg(W, Nit, cgform)

function [out, Elist]=pccg(W, Nit, cgform)
out=W;
alphat=9e-5;
g=getgrad(out);

for it=1:Nit
g_old = g;
g=getgrad(out);
if(it>1)
  printf('angle cosine : %f and cg test : %f and energy : %f \
	 \n',real(sum(diag(g'*d)))/(sqrt(real(sum(diag(d'*d)))*real(sum(diag(g'*g))))), \
	 real(sum(diag(g'*K(g_old))))/(sqrt(real(sum(diag(g'*g)))*real(sum(diag(K(g_old)'*K(g_old)))))), \
	 getE(out));
  if(cgform == 1)
    d=-K(g)+real(sum(diag(g'*K(g))))/(real(sum(diag(g_old'*K(g_old)))))*d;
  elseif(cgform == 2)
    d=-K(g)+real(sum(diag((g-g_old)'*K(g))))/(real(sum(diag((g-g_old)'*K(g_old)))))*d;
  elseif(cgform == 3)
    d=-K(g)+real(sum(diag((g-g_old)'*K(g))))/(real(sum(diag((g-g_old)'*d))))*d;
  end
else
  d=-K(g); 
end

gt=getgrad(out+alphat*d);
alpha=alphat*(real(sum(diag(g'*d))))/(real(sum(diag((g-gt)'*d))));
out=out+alpha*d;
Elist(it)=getE(out);
end

endfunction
