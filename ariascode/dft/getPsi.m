function [Psi, epsilon] = getPsi(W)

Y = W*sqrtm(inv(W'*O(W)));
mu = Y'*H(Y);

[U, epsilon] = eig(mu);
epsilon = real(diag(epsilon));

Psi = Y*U;

endfunction