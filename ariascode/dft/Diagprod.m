%# produces diag(a)*B more efficiently

function out = Diagprod(a,B)

out = (a*ones(1,size(B,2))).*B;

endfunction