function deriv = dx_dxi(H, xi, lambda, rho) 
  deriv = lambda * (lambda + (1 - lambda) * (H/xi)^rho)^(1/(rho) - 1);
end