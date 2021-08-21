function mom = calcmom(lambda, mu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta, A_0_tilde_no_delta_pz, ...
    A_1_tilde_no_delta_pz, c_0_tilde_no_delta_pz, ...
    c_1_tilde_no_delta_pz, p_z, calc_irfs, make_plots, A_1)

   % so I don't get a billion "singular" warnings 
   warning('off','all')

   shock_vars = A_1 * [xi_star; steady_state];
   
   xi_shock = shock_vars(1) + kappa;
   shock_state = shock_vars(2:end);

   H_inside = 0;
   
   agg_scale_factor = sqrt(n_periods / 5) * sqrt(omega * (1 - omega));
   irf_scale_factor = sqrt(n_periods / 5) * sqrt(omega * alpha / p_z * (1 - omega * alpha / p_z));
   
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
   
   y_irf = (log(Y_shock) - log(Y_star)) *  agg_scale_factor;
   
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
       
       s0_pz_trans = A_0_tilde_no_delta_pz(2:end, 2:end);
       s1_pz_trans = A_1_tilde_no_delta_pz(2:end, 2:end);
       s0_pz_intercept = c_0_tilde_no_delta_pz(2:end, :);
       s1_pz_intercept = c_1_tilde_no_delta_pz(2:end, :);

       
       [sswg, ss_dist_mat, ss_wage_mat] = calc_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage, low_wage, ...
                               wages_by_bin');
       
       % wage growth and absolute wage growth 
       % by bin, scaled by asymptotic variance of xi
       exposed_cfwg = calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_pz_trans, s1_pz_trans, ...
                               s0_pz_intercept, s1_pz_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage, low_wage, ...
                               wages_by_bin', ...
                               0);

       % this version shocks in the first period
       exposed_shockwg = calc_shock_wage_growth(s0_trans, ...
                                        s1_trans,s0_intercept,...
                                        s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                        s0_pz_intercept, s1_pz_intercept, omega,... 
                                        theta_grid, n_periods, ... 
                                        high_wage_tp1, low_wage_tp1, ...
                                        wages_by_bin', ...
                                        1);
       
       
       unexposed_cfwg = calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_pz_trans, s1_pz_trans, ...
                               s0_pz_intercept, s1_pz_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage, low_wage, ...
                               wages_by_bin', ...
                               0);

       % this version shocks in the first period
       % we calculate regular wage growth for the "shock" with 
       % the prevailing ending wages in in the post-shock period
       % because unexposed people can't fall down the ladder
       % last argument indicates fall probability in first 
       % period
       unexposed_shockwg = calc_shock_wage_growth(s0_trans, ...
                                        s1_trans,s0_intercept,...
                                        s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                        s0_pz_intercept, s1_pz_intercept, omega,... 
                                        theta_grid, n_periods, ... 
                                        high_wage_tp1, low_wage_tp1, ...
                                        wages_by_bin', ...
                                        0);
       
       quantile_dist_vec = reshape(ss_dist_mat * diag(steady_state), length(theta_grid)^2, 1);
       quantile_wage_vec = reshape(ss_wage_mat,length(theta_grid)^2, 1);
       targets = 0.01:0.01:0.99;
       q = zeros(size(targets));
       for i = 1:length(targets)
           q(i) = weighted_quantile(quantile_wage_vec, quantile_dist_vec, targets(i));
       end 
       
       % steady state wage growth percentile
       tenth_pctile = weighted_quantile(quantile_wage_vec, quantile_dist_vec, 0.1);
       
       exposed_cf_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage, low_wage, ...
                                  wages_by_bin', tenth_pctile, 0);

       exposed_shock_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_tp1, low_wage_tp1, ...
                                  wages_by_bin', tenth_pctile, 1);
                              
       unexposed_cf_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage, low_wage, ...
                                  wages_by_bin', tenth_pctile, 0);

       unexposed_shock_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_tp1, low_wage_tp1, ...
                                  wages_by_bin', tenth_pctile, 0);

       agg_shock_exposed_wg = zeros(5, 1);
       agg_noshock_exposed_wg = zeros(5, 1);
       agg_shock_unexposed_wg = zeros(5, 1);
       agg_noshock_unexposed_wg = zeros(5, 1);
       
       agg_shock_exposed_awg = zeros(5, 1);
       agg_noshock_exposed_awg = zeros(5, 1);
       agg_shock_unexposed_awg = zeros(5, 1);
       agg_noshock_unexposed_awg = zeros(5, 1);
       
       agg_shock_exposed_pctile = zeros(5, 1);
       agg_noshock_exposed_pctile = zeros(5, 1);
       agg_shock_unexposed_pctile = zeros(5, 1);
       agg_noshock_unexposed_pctile = zeros(5, 1);
       
       agg_prob = zeros(5, 1);
       
       % calculate wage growth and abs wage growth bins
       index = 1;
       cumulative_density = 0;

       while(cumulative_density < 0.25)
           cumulative_density = cumulative_density + steady_state(index);
           index = index + 1;
       end

       % correct last bit of loop
       index = index - 1;

       agg_shock_exposed_wg(1,1) = exposed_shockwg(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_exposed_wg(1,1) = exposed_cfwg(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));      
       agg_shock_unexposed_wg(1,1) = unexposed_shockwg(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_unexposed_wg(1,1) = unexposed_cfwg(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index))); 
       
       agg_shock_exposed_awg(1,1) = exposed_shockwg(1:index, 2)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_exposed_awg(1,1) = exposed_cfwg(1:index, 2)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));      
       agg_shock_unexposed_awg(1,1) = unexposed_shockwg(1:index, 2)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_unexposed_awg(1,1) = unexposed_cfwg(1:index, 2)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index))); 
       
       agg_shock_exposed_pctile(1,1) = exposed_shock_lt_pctile_probs(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_exposed_pctile(1,1) = exposed_cf_lt_pctile_probs(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));      
       agg_shock_unexposed_pctile(1,1) = unexposed_shock_lt_pctile_probs(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index)));
       agg_noshock_unexposed_pctile(1,1) = unexposed_cf_lt_pctile_probs(1:index, 1)' * ...
           (steady_state(1:index) ./ sum(steady_state(1:index))); 
       
       agg_prob(1, 1) = sum(steady_state(1:index));

       
       index = index + 1;
       next_start_index = index;
       while(cumulative_density < 0.5)
           cumulative_density = cumulative_density + steady_state(index);
           index = index + 1;
       end

       % correct last bit of loop
       index = index - 1;

       agg_shock_exposed_wg(2,1) = exposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_wg(2,1) = exposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_wg(2,1) = unexposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_wg(2,1) = unexposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_awg(2,1) = exposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_awg(2,1) = exposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_awg(2,1) = unexposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_awg(2,1) = unexposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_pctile(2,1) = exposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_pctile(2,1) = exposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_pctile(2,1) = unexposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_pctile(2,1) = unexposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 

       agg_prob(2, 1) = sum(steady_state(next_start_index:index));
       
       index = index + 1;
       next_start_index = index;


       while(cumulative_density < 0.75)
           cumulative_density = cumulative_density + steady_state(index);
           index = index + 1;
       end

       % correct last bit of loop
       index = index - 1;

       agg_shock_exposed_wg(3,1) = exposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_wg(3,1) = exposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_wg(3,1) = unexposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_wg(3,1) = unexposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_awg(3,1) = exposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_awg(3,1) = exposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_awg(3,1) = unexposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_awg(3,1) = unexposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_pctile(3,1) = exposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_pctile(3,1) = exposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_pctile(3,1) = unexposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_pctile(3,1) = unexposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 

       agg_prob(3, 1) = sum(steady_state(next_start_index:index));
       
       

       index = index + 1;
       next_start_index = index;

       while(cumulative_density < 0.95)
           cumulative_density = cumulative_density + steady_state(index);
           index = index + 1;
       end

       % correct last bit of loop
       index = index - 1;

       agg_shock_exposed_wg(4,1) = exposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_wg(4,1) = exposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_wg(4,1) = unexposed_shockwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_wg(4,1) = unexposed_cfwg(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_awg(4,1) = exposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_awg(4,1) = exposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_awg(4,1) = unexposed_shockwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_awg(4,1) = unexposed_cfwg(next_start_index:index, 2)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 
       
       agg_shock_exposed_pctile(4,1) = exposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_exposed_pctile(4,1) = exposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));      
       agg_shock_unexposed_pctile(4,1) = unexposed_shock_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index)));
       agg_noshock_unexposed_pctile(4,1) = unexposed_cf_lt_pctile_probs(next_start_index:index, 1)' * ...
           (steady_state(next_start_index:index) ./ sum(steady_state(next_start_index:index))); 

       agg_prob(4, 1) = sum(steady_state(next_start_index:index));
       
       
       index = index + 1;
       agg_shock_exposed_wg(5,1) = exposed_shockwg(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_exposed_wg(5,1) = exposed_cfwg(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));      
       agg_shock_unexposed_wg(5,1) = unexposed_shockwg(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_unexposed_wg(5,1) = unexposed_cfwg(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end))); 
       
       agg_shock_exposed_awg(5,1) = exposed_shockwg(index:end, 2)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_exposed_awg(5,1) = exposed_cfwg(index:end, 2)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));      
       agg_shock_unexposed_awg(5,1) = unexposed_shockwg(index:end, 2)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_unexposed_awg(5,1) = unexposed_cfwg(index:end, 2)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end))); 
       
       agg_shock_exposed_pctile(5,1) = exposed_shock_lt_pctile_probs(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_exposed_pctile(5,1) = exposed_cf_lt_pctile_probs(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));      
       agg_shock_unexposed_pctile(5,1) = unexposed_shock_lt_pctile_probs(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end)));
       agg_noshock_unexposed_pctile(5,1) = unexposed_cf_lt_pctile_probs(index:end, 1)' * ...
           (steady_state(index:end) ./ sum(steady_state(index:end))); 

       agg_prob(5, 1) = sum(steady_state(index:end));
       
       % which outcomes are associated with shocks
       S = [ones(5, 1); zeros(5, 1); ones(5, 1); zeros(5, 1);];
       % which outcomes are associated with shocks and exposure
       T = [zeros(5, 1); zeros(5, 1); ones(5, 1); ones(5, 1);].*S;
       % diagonal matrix of indicators for income bins and interactions of
       % income bins with T
       D = [repmat(diag(ones(5, 1)), 4, 1), repmat(diag(ones(5, 1)), 4, 1) .* T];
              
       y_vec_wg = [agg_shock_unexposed_wg; agg_noshock_unexposed_wg; ...
           agg_shock_exposed_wg; agg_noshock_exposed_wg];
       
       y_vec_awg = [agg_shock_unexposed_awg; agg_noshock_unexposed_awg; ...
           agg_shock_exposed_awg; agg_noshock_exposed_awg];
       
       y_vec_pctile = [agg_shock_unexposed_pctile; agg_noshock_unexposed_pctile; ...
           agg_shock_exposed_pctile; agg_noshock_exposed_pctile];
       
       e_xy_wg = zeros(11, 1);
       e_xy_awg = zeros(11, 1);
       e_xy_pctile = zeros(11, 1);
        
       probvec = [agg_prob; agg_prob; agg_prob; agg_prob];
       b =  alpha / p_z ;
       b_vec = b * T + (1 - b) * (1 - T);
       
       % compute the expectation E(y * S)
       % the only relevant state is S = 1 because the rest is 0
       e_xy_wg(1,:) = omega * (S .* b_vec .* probvec)' * y_vec_wg;
       e_xy_awg(1,:) =  omega * (S .* b_vec .* probvec)' * y_vec_awg;
       e_xy_pctile(1,:) =  omega * (S .* b_vec .* probvec)' * y_vec_pctile;
       
       % loop over the columns of the D matrix and
       % compute the expectation E(y * D)
       for i = 2:6
           D_i = D(:,i - 1);
           e_xy_wg(i,:) = agg_prob(i - 1) * (D_i .* (omega * (S .* b_vec) + ...
               (1 - omega) * ((1 - S) .* b_vec)))' * y_vec_wg;
           e_xy_awg(i,:) = agg_prob(i - 1) * (D_i .* (omega * (S .* b_vec) + ...
               (1 - omega) * ((1 - S) .* b_vec)))' * y_vec_awg;
           e_xy_pctile(i,:) = agg_prob(i - 1) * (D_i .* (omega * (S .* b_vec) + ...
               (1 - omega) * ((1 - S) .* b_vec)))' * y_vec_pctile;
       end
       
       % loop over the columns of the D matrix and
       % compute the expectation E(y * D * T)
       % T is already included in D
       for i = 7:11
           D_i = D(:,i - 1);
           e_xy_wg(i,:) = agg_prob(i - 6) * b * (D_i .* (omega * (S) + ...
               (1 - omega) * ((1 - S))))' * y_vec_wg;
           e_xy_awg(i,:) = agg_prob(i - 6) * b * (D_i .* (omega * (S) + ...
               (1 - omega) * ((1 - S))))' * y_vec_awg;
           e_xy_pctile(i,:) = agg_prob(i - 6) * b * (D_i .* (omega * (S) + ...
               (1 - omega) * ((1 - S))))' * y_vec_pctile;
       end
              
       state_stuff = [omega * agg_prob', omega * alpha / p_z * agg_prob']';
       
       xtx = [[omega; state_stuff], [state_stuff'; ...
           [[diag(agg_prob'); ...
           diag(omega * alpha / p_z * agg_prob')],...
           [diag(omega * alpha / p_z * agg_prob'); ...
           diag(omega * alpha / p_z * agg_prob')]]]];
       
       
       wage_growth = (xtx) \ e_xy_wg;
       wage_growth = wage_growth((end - 4):end) .* irf_scale_factor;
       abs_wage_growth = (xtx) \ e_xy_awg;
       abs_wage_growth = abs_wage_growth((end - 4):end) .* irf_scale_factor;
       
       tenth_pctile_probs = (xtx) \ e_xy_pctile;
       tenth_pctile_probs = tenth_pctile_probs((end - 4):end) .* irf_scale_factor;
       
   if make_plots > 0
       
      figure(6)
      
      
      subplot(2, 1, 1);
      bar([agg_shock_unexposed_wg, agg_shock_exposed_wg, agg_noshock_exposed_wg]);
      legend('Shock Unexposed','Shock Exposed','No Shock','Location','southoutside', ...
          'Orientation', 'horizontal');
      title("Breakdown Numbers")
      
      subplot(2, 1, 2);
      bar([agg_shock_unexposed_wg - agg_noshock_unexposed_wg,  ... 
          agg_shock_exposed_wg - agg_noshock_exposed_wg]);
      legend('Shock Unexposed - No Shock Unexposed',...
          'Shock Exposed - No Shock Exposed','Location','southoutside', ...
          'Orientation', 'horizontal');
      title("Diff Coefs")
      
      figure(7)
      subplot(2, 1, 1);
      bar([agg_shock_exposed_wg - agg_noshock_exposed_wg - ...
          (agg_shock_unexposed_wg - agg_noshock_unexposed_wg)]);
      title("Diff in Diff Coefs")
    
      subplot(2, 1, 2);
      bar(wage_growth);
      title("Regression Coefs")
   
   end
       
   else 
       wage_growth = zeros(5, 1);
       abs_wage_growth = zeros(5, 1);
       tenth_pctile_probs = zeros(5, 1);
   end
   
   lshare = (H_star * high_wage + L_star * low_wage) / Y_star;
   lshare_shock = (H_shock * high_wage_tp1 + L_shock * low_wage_tp1) / Y_shock;
      
   lshare_irf_sign = ((log(lshare_shock) - log(lshare)) > 0) * abs(log(lshare_shock) - log(lshare)) * 10;
   
   lshare_irf = (log(lshare_shock) - log(lshare)) * agg_scale_factor;
   
   premium = top_five_wages / bottom_five_wages;
   
   if make_plots > 0
       
      figure(5)
      theta_plot_x_axis = (high_wage * theta_grid + (low_wage) * (1 - theta_grid) ) ./ ...
          ((low_wage * (1 - min(theta_grid))) + high_wage * min(theta_grid));
      plot(cumsum(steady_state), theta_plot_x_axis, "-o");
      title("Density by Wage / W_l")
   end
   
   mom = [lshare; premium; y_irf; lshare_irf; y_irf_sign; lw_irf_sign; ...
       lshare_irf_sign; ...
       abs_wage_growth; ...
       wage_growth; ...
       tenth_pctile_probs];
end