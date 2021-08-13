
function Y = calc_Y(H, L, X, mu, sigma, v, H_inside) 
  if H_inside == 0
    Y = ((mu * (H^sigma) + (1 - mu) * (X ^ sigma))^(1/sigma))^v;
  else
    Y = ((mu * (L^sigma) + (1 - mu) * (X ^ sigma))^(1/sigma))^v;
  end
end
