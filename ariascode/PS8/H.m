%# calculates the LHS of the equation (in this case Schrodinger eqn)

function out = H(W)

global gbl_Vdual;

out = -0.5*L(W) + cIdag(Diagprod(gbl_Vdual,cI(W)));

endfunction