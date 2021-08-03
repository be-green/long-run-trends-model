
function Y = calc_Y(H, X, mu, sigma, v) 
  Y = ((mu * (H^sigma) + (1 - mu) * (X ^ sigma))^(1/sigma))^v;
end
