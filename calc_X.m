
function X = calc_X(xi, L, lambda, rho) 
  X = (lambda * (xi^rho) + (1 - lambda) * (L ^ rho)) ^ (1 / rho);
end