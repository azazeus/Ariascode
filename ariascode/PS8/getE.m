function out = getE(W)

global gbl_Vdual;

N=W'*O(W);
A=cI(W)*inv(N);
B=cI(W);
n=diagouter(A,B);

out=real(-0.5*sum(diag(W'*L(W*inv(N)))) + gbl_Vdual'*n);

endfunction