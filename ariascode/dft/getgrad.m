%# calculates gradients

function grad=getgrad(W)

global gbl_f;

N=W'*O(W);
grad = gbl_f*(H(W) - (O(W*inv(N)))*(W'*H(W)))*inv(N);

endfunction