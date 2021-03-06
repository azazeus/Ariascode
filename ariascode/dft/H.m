%# calculates the LHS of the equation (in this case Schrodinger eqn)

function out = H(W)

global gbl_Vdual;
global gbl_f;

N=W'*O(W);
A=cI(W)*inv(N);
B=cI(W);
n=gbl_f*diagouter(A,B);
out = -0.5*L(W) + \
    cIdag(Diagprod((gbl_Vdual+cJdag(O(Linv(-4*pi*O(cJ(n))))) + \
		    cJdag(O(cJ(excVWN(n))))+ Diagprod(excpVWN(n), cJdag(O(cJ(n))))),cI(W)));

endfunction