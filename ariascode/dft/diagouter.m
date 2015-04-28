%# calculates diagonal elements of A*B_dag

function out= diagouter(A,B)

out=sum(A.*conj(B),2);

endfunction