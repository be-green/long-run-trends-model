function deriv = dx_dh (H, xi, lambda,  rho) 
  deriv = (1 - lambda) / lambda * (H/xi)^(rho - 1) * dx_dxi(H, xi, lambda, rho);
end