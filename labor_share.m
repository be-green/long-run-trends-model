function lshare = labor_share(L, X, H, xi, Y, lambda, mu, rho, sigma) 
  lshare = 1 - (lambda) * (1 - mu) * (mu * (L/X)^sigma + (1 - mu))^((1 / sigma) - 1) * ...
  (lambda + (1 - lambda) * (H/xi)^rho)^((1/rho) - 1) * xi / Y;
end