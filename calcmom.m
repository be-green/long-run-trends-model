function mom = calcmom(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta, calc_irfs)

   xi_shock = xi_star + kappa;
   shock_state = steady_state - alpha * steady_state;
   shock_state(1) = shock_state(1) + alpha;

   H_inside = 0;
   
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
   X_star = calc_X(xi_star, H_star, L_star, lambda, rho, H_inside);
   Y_star = calc_Y(H_star, L_star, X_star, mu, sigma, v, H_inside);
   
   % high and low wages at steady state
   high_wage = w_h(H_star, L_star, xi_star, rho, sigma, mu, lambda, v, H_inside);
   low_wage = w_l(H_star, L_star, xi_star, rho, sigma, mu, lambda, v, H_inside);
   
   % shock state production values
   X_shock = calc_X(xi_shock, H_shock, L_shock, lambda, rho, H_inside);
   Y_shock = calc_Y(H_shock, L_shock, X_shock, mu, sigma, v, H_inside);
   
   % high and low wages at shock state
   high_wage_tp1 = w_h(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v, H_inside);
   low_wage_tp1 = w_l(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v, H_inside);
   
   % signs of IRF for optimizer
   y_irf_sign = ((Y_shock - Y_star) < 0) * abs(Y_shock - Y_star) * 10;
   lw_irf_sign = ((low_wage_tp1 - low_wage) > 0) * abs(low_wage_tp1 - low_wage) * 10;
   
   y_irf = (log(Y_shock) - log(Y_star)) /  (kappa / sqrt(xi_var));
   
   % wages by each part of theta grid
   % for steady state ratio
   wages_by_bin = high_wage .* theta_grid + low_wage .* (1 - theta_grid);   
   
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
   
   if calc_irfs > 0
   
       s0_trans = A_0_tilde_no_delta(2:end, 2:end);
       s1_trans = A_1_tilde_no_delta(2:end, 2:end);
       s0_intercept = c_0_tilde_no_delta(2:end, :);
       s1_intercept = c_1_tilde_no_delta(2:end, :);

       % wage growth and absolute wage growth 
       % by bin, scaled by asymptotic variance of xi
       [sswg, ss_dist_mat, ss_wage_mat] = calc_wage_growth(s0_trans, ...
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
       ss_wage_growth = sswg(:,1)/ (kappa / sqrt(xi_var));
       ss_abs_wage_growth = sswg(:,2)/ (kappa / sqrt(xi_var));
%        figure(3)
%        bar(log(theta_grid), ss_abs_wage_growth)
       
       % wage growth following shock by income bin
       shock_wage_growth = shockwg(:,1);
       shock_abs_wage_growth = shockwg(:,2);

       shock_wage_growth = shock_wage_growth /  (kappa / sqrt(xi_var));
       shock_abs_wage_growth = shock_abs_wage_growth /  (kappa / sqrt(xi_var));

       % marginal effects by bin
       % since this is a linear combo I think this is the same
       % as aggregating and then calculating
       wage_growth = shock_wage_growth - ss_wage_growth;
       abs_wage_growth = shock_abs_wage_growth - ss_abs_wage_growth;
       
       quantile_dist_vec = reshape(ss_dist_mat * diag(steady_state), length(theta_grid)^2, 1);
       quantile_wage_vec = reshape(ss_wage_mat,length(theta_grid)^2, 1);
       targets = 0.01:0.01:0.99;
       q = zeros(size(targets));
       for i = 1:length(targets)
           q(i) = weighted_quantile(quantile_wage_vec, quantile_dist_vec, targets(i));
       end 
       
       tenth_pctile = weighted_quantile(quantile_wage_vec, quantile_dist_vec, 0.25);
       ss_lt_pctile_probs = calc_wage_growth_lt_pctile(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage, low_wage_tp1, ...
                                  wages_by_bin', tenth_pctile);

       shock_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_tp1, low_wage_tp1, ...
                                  wages_by_bin', tenth_pctile);

       wage_growth = wage_growth';
       abs_wage_growth = abs_wage_growth';

    %    wage_growth = shock_wage_growth - ss_wage_growth;
       tenth_pctile_effects = (shock_lt_pctile_probs - ss_lt_pctile_probs) /  (kappa / sqrt(xi_var));
       tenth_pctile_effects = tenth_pctile_effects';


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

       bottom_twentyfive_tenthpctile = tenth_pctile_effects(1:index) * ...
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

       twentyfive_to_fifty_tenthpctile = tenth_pctile_effects(next_start_index:index) * ...
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

      fifty_to_seventyfive_tenthpctile = tenth_pctile_effects(next_start_index:index) * ...
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

       seventyfive_to_ninetyfive_tenthpctile = tenth_pctile_effects(next_start_index:index) * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));

       index = index + 1;

      top_five_wagegrowth = wage_growth(index:end) * ...
          (steady_state(index:end) ./ sum(steady_state(index:end)));

      top_five_abswagegrowth = abs_wage_growth(index:end) * ...
       (steady_state(index:end) ./ sum(steady_state(index:end)));

      top_five_tenthpctile = tenth_pctile_effects(index:end) * ...
       (steady_state(index:end) ./ sum(steady_state(index:end)));

        wage_growth = [bottom_twentyfive_wagegrowth; twentyfive_to_fifty_wagegrowth;...
           fifty_to_seventyfive_wagegrowth; seventyfive_to_ninetyfive_wagegrowth; ...
           top_five_wagegrowth];

       abs_wage_growth = [bottom_twentyfive_abswagegrowth; twentyfive_to_fifty_abswagegrowth;...
           fifty_to_seventyfive_abswagegrowth; seventyfive_to_ninetyfive_abswagegrowth; ...
           top_five_abswagegrowth];

       tenth_pctile_probs = [bottom_twentyfive_tenthpctile; twentyfive_to_fifty_tenthpctile;...
           fifty_to_seventyfive_tenthpctile; seventyfive_to_ninetyfive_tenthpctile; ...
           top_five_tenthpctile];
   else 
       wage_growth = zeros(5, 1);
       abs_wage_growth = zeros(5, 1);
       tenth_pctile_probs = zeros(5, 1);
   end
   
   lshare = (H_star * high_wage + L_star * low_wage) / Y_star;
   lshare_shock = (H_shock * high_wage_tp1 + L_shock * low_wage_tp1) / Y_shock;
      
   lshare_irf_sign = ((log(lshare_shock) - log(lshare)) > 0) * abs(log(lshare_shock) - log(lshare)) * 10;
   
   lshare_irf = (log(lshare_shock) - log(lshare)) / (kappa / sqrt(xi_var));
   
   premium = top_five_wages / bottom_five_wages;
   mom = [lshare; premium; y_irf; lshare_irf; y_irf_sign; lw_irf_sign; ...
       lshare_irf_sign; ...
       abs_wage_growth; ...
       wage_growth; ...
       tenth_pctile_probs];
end