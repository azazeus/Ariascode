% Performs finite difference test of getE() and getgrad()
% 
% Usage: fdtest(W,S)
%
% W: starting point for test (size: prod(S) x Ns)
% S: Dimensions of 3d data

function fdtest(W)
  %# Compute intial energy and gradient
  E0=getE(W)
  g0=getgrad(W);

  %# Choose a random direction to explore
  dW=randn(size(W))+i*randn(size(W));

  %# Explore a range of step sizes decreasing by powers of ten
  for delta=10.^[1:-1:-9]
    %# Directional derivative formula
    dE=2*real(trace(g0'*delta*dW));
    
    %# Print ratio of actual change to expected change, along with
     # estimate
    %#   of the error in this quantity due to rounding
    printf(' %20.16f\n %20.16f\n\n', ...
    (getE(W+delta*dW)-E0)/dE, sqrt(size(W,1))*eps/abs(dE) );
  end
endfunction