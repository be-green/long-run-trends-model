
function X = calc_X(xi, H, L, lambda, rho, H_inside) 
  if H_inside == 0
    X = (lambda * (xi^rho) + (1 - lambda) * (L ^ rho)) ^ (1 / rho);
  else
    X = (lambda * (xi^rho) + (1 - lambda) * (H ^ rho)) ^ (1 / rho);  
  end 
end