function mom = calcmom(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v)

   xi_shock = xi_star + kappa;
   shock_state = steady_state - alpha * steady_state;
   shock_state(1) = shock_state(1) + alpha;

   H_star = theta_grid * steady_state;
   L_star = 1 - H_star;
   
   shock_vec = [xi_shock; shock_state(1:(end - 1))];
   for i=1:(n_periods - 1)
       shock_vec = (1 - omega) * (A_0_tilde * shock_vec + c_0_tilde) + ...
           omega  * (A_1_tilde * shock_vec + c_1_tilde);
   end
   
   shock_vec = [shock_vec; 1 - sum(shock_vec(2:end))];
   
   xi_shock = shock_vec(1,:);
   shock_state = shock_vec(2:end,:);
   
   H_shock = theta_grid * shock_state;
   L_shock = 1 - H_shock;
   
   % translates to 0-1 constraint from unconstrained space
   lambda = normcdf(x(1, 1));
   mu = normcdf(x(2, 1));
   
   % steady state production values
   X_star = calc_X(xi_star, L_star, lambda, rho);
   Y_star = calc_Y(H_star, X_star, mu, sigma, v);
   
   % high and low wages at steady state
   high_wage = w_h(H_star, L_star, xi_star, rho, sigma, mu, lambda, v);
   low_wage = w_l(H_star, L_star, xi_star, rho, sigma, mu, lambda, v);
   
   % shock state production values
   X_shock = calc_X(xi_shock, L_shock, lambda, rho);
   Y_shock = calc_Y(H_shock, X_shock, mu, sigma, v);
   
   % high and low wages at shock state
   high_wage_tp1 = w_h(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v);
   low_wage_tp1 = w_l(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v);
   
   % signs of IRF for optimizer
   y_irf_sign = sign(Y_shock - Y_star);
   lw_irf_sign = sign(low_wage_tp1 - low_wage);
   
   y_irf = (log(Y_shock) - log(Y_star)) / sqrt(xi_var) * kappa;
   
   % wages by each part of theta grid
   % for steady state ratio
   wages_by_bin = high_wage .* theta_grid + low_wage .* (1 - theta_grid);   
   
   s0_trans = A_0_tilde(2:end, 2:end);
   s1_trans = A_1_tilde(2:end, 2:end);
   s0_intercept = c_0_tilde(2:end, :);
   s1_intercept = c_1_tilde(2:end, :);
   
   % wage growth and absolute wage growth 
   % by bin, scaled by asymptotic variance of xi
   sswg = calc_wage_growth(s0_trans, ...
                           s1_trans, s0_intercept, ...
                           s1_intercept, omega, ...
                           theta_grid, n_periods, ...
                           high_wage, low_wage, ...
                           wages_by_bin');
                                                 
   % this version shocks in the first period
   shockwg = calc_shock_wage_growth(s0_trans, ...
                                    s1_trans,s0_intercept,...
                                    s1_intercept, omega,... 
                                    theta_grid, n_periods, ... 
                                    high_wage_tp1, low_wage_tp1, ...
                                    wages_by_bin');

   % average wage growth and abs wage growth at ergodic steady state
   ss_wage_growth = sswg(:,1)/sqrt(xi_var) * kappa;
   ss_abs_wage_growth = sswg(:,2)/sqrt(xi_var) * kappa;
   
   % wage growth following shock by income bin
   shock_wage_growth = shockwg(:,1);
   shock_abs_wage_growth = shockwg(:,2);
   
   shock_wage_growth = shock_wage_growth / sqrt(xi_var) * kappa;
   shock_abs_wage_growth = shock_abs_wage_growth / sqrt(xi_var) * kappa;
   
   % marginal effects by bin
   % since this is a linear combo I think this is the same
   % as aggregating and then calculating
   wage_growth = shock_wage_growth - ss_wage_growth;
   abs_wage_growth = shock_abs_wage_growth - ss_abs_wage_growth;
   
   wage_growth = wage_growth';
   abs_wage_growth = abs_wage_growth';
   
   % calculate wage premium at steady state
   bottom_five = 1;
   cumulative_density = 0;
   
   while(cumulative_density < 0.25)
       cumulative_density = cumulative_density + steady_state(bottom_five);
       bottom_five = bottom_five + 1;
   end
   
   bottom_five_wages = wages_by_bin(bottom_five - 1);
   
   top_five = size(steady_state,1);
   cumulative_density = 0;
   
   while(cumulative_density < 0.25)
       cumulative_density = cumulative_density + steady_state(top_five);
       top_five = top_five - 1;
   end
   
   top_five_wages = wages_by_bin(top_five + 1);
   
   % calculate wage growth and abs wage growth bins
   index = 1;
   cumulative_density = 0;
   
   while(cumulative_density < 0.25)
       cumulative_density = cumulative_density + steady_state(index);
       index = index + 1;
   end
   
   % correct last bit of loop
   index = index - 1;
   
   bottom_twentyfive_wagegrowth = wage_growth(1:index) * ...
       (steady_state(1:index) ./ sum(steady_state(1:index)));
   
   bottom_twentyfive_abswagegrowth = abs_wage_growth(1:index) * ...
       (steady_state(1:index) ./ sum(steady_state(1:index)));
   
   index = index + 1;
   next_start_index = index;
   while(cumulative_density < 0.5)
       cumulative_density = cumulative_density + steady_state(index);
       index = index + 1;
   end
   
   % correct last bit of loop
   index = index - 1;
   
   twentyfive_to_fifty_wagegrowth = wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   twentyfive_to_fifty_abswagegrowth = abs_wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   index = index + 1;
   next_start_index = index;
   
   
   while(cumulative_density < 0.75)
       cumulative_density = cumulative_density + steady_state(index);
       index = index + 1;
   end
   
   % correct last bit of loop
   index = index - 1;
   
   fifty_to_seventyfive_wagegrowth = wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   fifty_to_seventyfive_abswagegrowth = abs_wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   index = index + 1;
   next_start_index = index;
   
   while(cumulative_density < 0.95)
       cumulative_density = cumulative_density + steady_state(index);
       index = index + 1;
   end
   
   % correct last bit of loop
   index = index - 1;
   
   seventyfive_to_ninetyfive_wagegrowth = wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   seventyfive_to_ninetyfive_abswagegrowth = abs_wage_growth(next_start_index:index) * ...
       (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
   
   index = index + 1;
   
  top_five_wagegrowth = wage_growth(index:end) * ...
      (steady_state(index:end) ./ sum(steady_state(index:end)));
   
  top_five_abswagegrowth = abs_wage_growth(index:end) * ...
   (steady_state(index:end) ./ sum(steady_state(index:end)));

    wage_growth = [bottom_twentyfive_wagegrowth; twentyfive_to_fifty_wagegrowth;...
       fifty_to_seventyfive_wagegrowth; seventyfive_to_ninetyfive_wagegrowth; ...
       top_five_wagegrowth];

   abs_wage_growth = [bottom_twentyfive_abswagegrowth; twentyfive_to_fifty_abswagegrowth;...
       fifty_to_seventyfive_abswagegrowth; seventyfive_to_ninetyfive_abswagegrowth; ...
       top_five_abswagegrowth];
   
   lshare = (H_star * high_wage + L_star * low_wage) / Y_star;
   lshare_shock = (H_shock * high_wage_tp1 + L_shock * low_wage_tp1) / Y_shock;
   
   lshare_irf_sign = sign(log(lshare_shock) - log(lshare));
   
   lshare_irf = (log(lshare_shock) - log(lshare)) / sqrt(xi_var) * kappa;
   
   premium = top_five_wages / bottom_five_wages;
   mom = [lshare; premium; y_irf; lshare_irf; y_irf_sign; lw_irf_sign; ...
       lshare_irf_sign; ...
       abs_wage_growth; ...
       wage_growth];
end