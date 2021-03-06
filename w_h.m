
function high_wage = w_h(H, L, xi, rho, sigma, mu, lambda, v, H_inside)
  if H_inside == 0
      high_wage = ((1 - mu) * (lambda * xi ^ rho - ...
          L ^ rho * (-1 + lambda)) ^ (sigma / rho) + ...
          mu * H ^ sigma) ^ (v / sigma) * v * ...
          mu * H ^ sigma / H / ((1 - mu) * ...
          (lambda * xi ^ rho - L ^ rho * ...
          (-1 + lambda)) ^ (sigma / rho) + ...
          mu * H ^ sigma);
  else
      high_wage = ((1 - mu) * ...
          (lambda * xi ^ rho - H ^ rho * (-1 + lambda)) ...
          ^ (sigma / rho) + ...
          mu * L ^ sigma) ^ (v / sigma) * ...
          v * (-1 + lambda) * (lambda * xi ^ rho - H ^ rho * ...
          (-1 + lambda)) ^ (sigma / rho) * (-1 + mu) * H ^ rho / ...
          (lambda * xi ^ rho + (1 - lambda) * H ^ rho) / H / ((1 - mu) * ...
          (lambda * xi ^ rho - H ^ rho * (-1 + lambda)) ^ (sigma / rho) + ...
          mu * L ^ sigma);
  end
end


