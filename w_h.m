
function high_wage = w_h(H, L, xi, rho, sigma, mu, lambda, v)
  high_wage = ((1 - mu) * (lambda * xi ^ rho - L ^ rho * (-1 + lambda)) ^ (sigma / rho) + mu * H ^ sigma) ^ (v / sigma) * v * mu * H ^ sigma / H / ((1 - mu) * (lambda * xi ^ rho - L ^ rho * (-1 + lambda)) ^ (sigma / rho) + mu * H ^ sigma);
end

