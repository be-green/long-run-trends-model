function [mom, diff_coefs, reg_coefs] = calcmom_no_displacement(lambda, mu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0, A_1, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_no_delta, A_1_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta, A_0_no_delta_pz, ...
    A_1_no_delta_pz, c_0_tilde_no_delta_pz, ...
    c_1_tilde_no_delta_pz, p_z, calc_irfs, make_plots, ...
    scale_period, H_inside)

   % scaling factors to convert from one shock units to 1 SD units
   agg_scale_factor = sqrt(scale_period) * sqrt(omega * (1 - omega));
   irf_scale_factor = sqrt(scale_period) * sqrt(omega * alpha / p_z * (1 - omega * alpha / p_z));
   % steady state levels of H and L
   H_star = theta_grid * steady_state;
   L_star = 1 - H_star;
   
   % calculate steady state aggregates & wages
   
   % steady state production values
   X_star = calc_X(xi_star, H_star, L_star, lambda, rho, H_inside);
   Y_star = calc_Y(H_star, L_star, X_star, mu, sigma, v, H_inside);
   
   % high and low wages at steady state
   high_wage_ss = w_h(H_star, L_star, xi_star, rho, sigma, mu, lambda, v, H_inside);
   low_wage_ss = w_l(H_star, L_star, xi_star, rho, sigma, mu, lambda, v, H_inside);
   
   % xi has an intercept that we need to include
   % that isn't just in the A_0 and A_1 matrices
   % grab that from the intercept terms in the VECM representation
   % since they are the same in both cases.
   int_for_xi_0 = zeros(size(A_0, 1), 1);
   int_for_xi_0(1) = c_0_tilde(1);
   int_for_xi_1 = zeros(size(A_0, 1), 1);
   int_for_xi_1(1) = c_1_tilde(1);
   
   A_1_no_displacement = A_0;
   
   % get single shock distribution for wage calculations
   shock_vars = A_1_no_displacement * [xi_star; steady_state] + int_for_xi_1;
   
   for i=1:(n_periods - 1)
       shock_vars = (1 - omega) * (A_0 * shock_vars + int_for_xi_0) + ...
           omega  * (A_1 * shock_vars + int_for_xi_1);
   end
   
   %Note: if we add a nonzero mass with zero skill, almost nothing would
   %change, except that we would have 1-p0 in place of 1 here. 
   %shock_vec = [shock_vec; 1 - sum(shock_vec(2:end))];
   
   % We also might want to put the zero in on the grid somewhere, in which case we
   %need to tweak things to be the following
   % shock_vec = [shock_vec(1); p0; shock_vec(2:end); 1-p0 - sum(shock_vec(2:end))];
   xi_shock = shock_vars(1,:);
   shock_state = shock_vars(2:end,:);
   
   % H & L after a shock + n_periods - 1 iteration
   H_shock = theta_grid * shock_state; %if we add the zero to theta_grid, this doesn't change
   L_shock = 1 - H_shock;
  
   % shock state production values
   X_shock = calc_X(xi_shock, H_shock, L_shock, lambda, rho, H_inside);
   Y_shock = calc_Y(H_shock, L_shock, X_shock, mu, sigma, v, H_inside);
   
   % high and low wages at shock state
   high_wage_shock = w_h(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v, H_inside);
   low_wage_shock = w_l(H_shock, L_shock, xi_shock, rho, sigma, mu, lambda, v, H_inside);
   
   % signs of IRF for optimizer
   y_irf_sign = ((Y_shock - Y_star) < 0) * abs(Y_shock - Y_star) * 10;
   lw_irf_sign = ((low_wage_shock - low_wage_ss) > 0) * abs(low_wage_shock - low_wage_ss) * 10;
   
   y_irf = (log(Y_shock) - log(Y_star)) *  agg_scale_factor;
   
   % Do the same thing as abovefor no shock 
   % single no shock period
   no_shock_vars = A_0 * [xi_star; steady_state];
   
   % add back intercept
   xi_no_shock = no_shock_vars(1) + c_0_tilde(1);
   no_shock_state = no_shock_vars(2:end);
   
   % remove top of grid for iteration
   no_shock_vec = [xi_no_shock; no_shock_state(1:(end))];
   
   % iterate forward the rest of n_periods
   for i=1:(n_periods - 1)
       no_shock_vec = (1 - omega) * (A_0 * no_shock_vec + int_for_xi_0) + ...
           omega  * (A_1 * no_shock_vec + int_for_xi_1);
   end
   % no_shock_vec = [no_shock_vec; 1 - sum(no_shock_vec(2:end))];

   xi_no_shock = no_shock_vec(1,:);
   no_shock_state = no_shock_vec(2:end,:);
   
   
   % H & L after a shock + n_periods - 1 iteration
   H_no_shock = theta_grid * no_shock_state; %if we add the zero to theta_grid, this doesn't change
   L_no_shock = 1 - H_no_shock;
   
   % no shock state production values
   % not used for now but you never know
   X_no_shock = calc_X(xi_no_shock, H_no_shock, L_no_shock, lambda, rho, H_inside);
   Y_no_shock = calc_Y(H_no_shock, L_no_shock, X_no_shock, mu, sigma, v, H_inside);
   
   % high and low wages at no shock state
   high_wage_no_shock = w_h(H_no_shock, L_no_shock, xi_no_shock, rho, sigma, mu, lambda, v, H_inside);
   low_wage_no_shock = w_l(H_no_shock, L_no_shock, xi_no_shock, rho, sigma, mu, lambda, v, H_inside);
   
   % display some stuff because it's helpful
   if(make_plots>0)
      disp(['High Wage Post Shock: ', num2str(high_wage_shock)])
      disp(['Low Wage Post Shock: ', num2str(low_wage_shock)])
      disp(['High Wage Post No Shock: ', num2str(high_wage_no_shock)])
      disp(['Low Wage Post No Shock: ', num2str(low_wage_no_shock)])
   end
   
   % wages by each part of theta grid
   % for steady state ratio
   wages_by_bin = high_wage_ss .* theta_grid + low_wage_ss .* (1 - theta_grid);   
   
   % calculate wage premium at steady state
   % this is defined as 75% wages / 25% wages
   bottom_twentyfive = 1;
   cumulative_density = 0;
   
   % iterate up until reaching point above 25% of mass
   while(cumulative_density < 0.25)
       cumulative_density = cumulative_density + steady_state(bottom_twentyfive);
       bottom_twentyfive = bottom_twentyfive + 1;
   end
   
   % go back a single cell and calc wages
   % max operator is in case of mass points at bottom
   bottom_twentyfive_wages = wages_by_bin(max(bottom_twentyfive - 1,1));
   
   % same idea but iterate backwards from the top
   % until you reach the 75%
   top_twentyfive = size(steady_state,1);
   cumulative_density = 0;
   
   while(cumulative_density < 0.25)
       cumulative_density = cumulative_density + steady_state(top_twentyfive);
       top_twentyfive = top_twentyfive - 1;
   end
   
   % min operator is in case of mass points
   top_twentyfive_wages = wages_by_bin(min(top_twentyfive + 1,end));
   
   if calc_irfs > 0
   
       % the transition matrices for the theta grid have no 
       % death, so we grab all values except the xi iterations
       % for the purposes of shocks
       s0_trans = A_0_no_delta(2:end, 2:end);
       s1_trans = A_1_no_delta(2:end, 2:end);
       s0_intercept = int_for_xi_0(2:end, :);
       s1_intercept = int_for_xi_1(2:end, :);
       
       % this is the same thing as above except w/ the fall probability
       % set to pz instead of alpha
       s0_pz_trans = A_0_no_delta_pz(2:end, 2:end);
       s1_pz_trans = A_1_no_delta_pz(2:end, 2:end);
       s0_pz_intercept = int_for_xi_0(2:end, :);
       s1_pz_intercept = int_for_xi_0(2:end, :);

       
       % calculate wage growth given a shock in the first period
       [sswg_shockwage, ss_dist_mat_shockwage, ss_wage_mat_shockwage] = ...
        calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_trans, s1_trans, ...
                               s0_intercept, s1_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage_shock, low_wage_shock, ...
                               wages_by_bin', ...
                               0); % no displacement first period
       
       % calculate wage growth given no shock in the first period
       [sswg_noshockwage, ss_dist_mat_noshockwage, ss_wage_mat_noshockwage] =  ...
        calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_trans, s1_trans, ...
                               s0_intercept, s1_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage_no_shock, low_wage_no_shock, ...
                               wages_by_bin', ...
                               0); % no shock first period
                           
       % steady state wage growth & abs_wage growth, 
       % taking into account wage changes given
       % shocks
       sswg = omega * sswg_shockwage + (1 - omega) * sswg_noshockwage;
       
       % steady state distributions and wages used for the 10th percentil
       % calculations later
       ss_dist_mat = omega * ss_dist_mat_shockwage + (1 - omega) * ss_dist_mat_noshockwage;
       ss_wage_mat = omega * ss_wage_mat_shockwage + (1 - omega) * ss_wage_mat_noshockwage;
       
       % expected wage growth & abs wage growth from the steady state
       expected_wage_growth = steady_state' * sswg(:,1);
       expected_abs_wage_growth = steady_state' * sswg(:,2);
                           
       % wage growth and absolute wage growth 
       
       % calculates the counterfactual for the exposed group's wage growth
       % and absolute wage growth. omega = 0 in the first period, meaning
       % there is no shock. after that, omega is normal
       exposed_cfwg = calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_pz_trans, s1_pz_trans, ...
                               s0_pz_intercept, s1_pz_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage_no_shock, low_wage_no_shock, ...
                               wages_by_bin', ...
                               0); % shock

       % this version shocks in the first period, uses the exposure fall
       % probabilites for that period
       exposed_shockwg = calc_shock_wage_growth(s0_trans, ...
                                        s1_trans,s0_intercept,...
                                        s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                        s0_pz_intercept, s1_pz_intercept, omega,... 
                                        theta_grid, n_periods, ... 
                                        high_wage_shock, low_wage_shock, ...
                                        wages_by_bin', ...
                                        0);
       
       % this version does not shock in the first period
       unexposed_cfwg = calc_shock_wage_growth(s0_trans, ...
                               s1_trans, s0_intercept, ...
                               s1_intercept,s0_pz_trans, s1_pz_trans, ...
                               s0_pz_intercept, s1_pz_intercept, omega, ...
                               theta_grid, n_periods, ...
                               high_wage_no_shock, low_wage_no_shock, ...
                               wages_by_bin', ...
                               0);

       % this version looks like there is "no shock" in the first period,
       % but that's just because there is no fall probability. It is using
       % the prevailing wages after a shock, indicating that in aggregate
       % the human capital has moved.
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
                                        high_wage_shock, low_wage_shock, ...
                                        wages_by_bin', ...
                                        0);
       
       % grab the aggregate steady state wage response matrix
       % weight it by the steady state density at each theta gridpoint
       % then re-shape it to calculate a single vector to calculate a
       % weighted quantile
       ss_dist_vec = reshape(ss_dist_mat * diag(steady_state), length(theta_grid)^2, 1);
       shock_dist_vec = reshape(ss_dist_mat_shockwage * diag(steady_state), length(theta_grid)^2, 1);
       noshock_dist_vec = reshape(ss_dist_mat_noshockwage * diag(steady_state), length(theta_grid)^2, 1);
       
       % this calculates E(WG) and E(WG^2) for every eventuality
       % the first of these is also used in calculating the quantiles of
       % wages
       % the second used in calculating variances via E(X^2) - E(X)^2
       % we first calculate X^2 - X for each paths in our grid
       % and then weight them by the probability associated w/ each
       % gridpoint, calculated as the P(starting_point) * P(ending_point |
       % starting_point) and subsequently aggregating
       first_moments = reshape(ss_wage_mat,length(theta_grid)^2, 1);
       second_moments = reshape(ss_wage_mat.^2,length(theta_grid)^2, 1);
       aggregate_variance = second_moments' * ss_dist_vec  - (first_moments' * ss_dist_vec)^2;
       
       first_moments_shock = reshape(ss_wage_mat_noshockwage,length(theta_grid)^2, 1);
       first_moments_noshock = reshape(ss_wage_mat_shockwage,length(theta_grid)^2, 1);
       
       % sd is sqrt var
       aggregate_sd = sqrt(aggregate_variance);
       
       targets = 0.01:0.01:0.99;
       q = zeros(size(targets));
       for i = 1:length(targets)
           q(i) = weighted_quantile(first_moments, ss_dist_vec, targets(i));
       end 
       
       if make_plots > 0
           figure
           plot(targets,q)
           title('Quantiles of wage growth distribution');
           xline(0.05);
           xline(0.25);
           xline(0.5);
           xline(0.75);
           xline(0.95);
           
           %TODO: we could add lines by income group. This might be a
           %pretty useful diagnostic to make sure things don't work for the
           %wrong reason
           
       end
       
       
       % steady state wage growth percentile
       shock_tenth_pctile = weighted_quantile(first_moments_shock, shock_dist_vec, 0.1);
       noshock_tenth_pctile = weighted_quantile(first_moments_noshock, noshock_dist_vec, 0.1);

       % this calculates the probability of falling below the 10th
       % percentile given that you are exposed and there is no shock
       exposed_cf_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_no_shock, low_wage_no_shock, ...
                                  wages_by_bin', noshock_tenth_pctile, 0);

       % this calculates the probability of falling below the 10th
       % percentile given that you are exposed and there is a shock
       exposed_shock_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_shock, low_wage_shock, ...
                                  wages_by_bin', shock_tenth_pctile, 0);
       
       % this calculates the probability of falling below the 10th
       % percentile given that you are unexposed and there is no shock                              
       unexposed_cf_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_no_shock, low_wage_no_shock, ...
                                  wages_by_bin', noshock_tenth_pctile, 0);
                              
       % this calculates the probability of falling below the 10th
       % percentile given that you are unexposed and there is a shock
       unexposed_shock_lt_pctile_probs = calc_wage_growth_lt_pctile_shock(s0_trans, ...
                                  s1_trans,s0_intercept,...
                                  s1_intercept, s0_pz_trans, s1_pz_trans, ...
                                  s0_pz_intercept, s1_pz_intercept, omega,... 
                                  theta_grid, n_periods, ... 
                                  high_wage_shock, low_wage_shock, ...
                                  wages_by_bin', shock_tenth_pctile, 0);

                              
       % pre allocated for aggregated versions of the responses
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
       quantile_targets = [0.25,0.5,0.75,0.95,1];
       
       ss_wg_by_income = zeros(5, 1);
       ss_awg_by_income = zeros(5, 1);
       
       % calculate wage growth and abs wage growth bins
       % indexing requires the vector of CDF levels
       ss_cdf = [cumsum(steady_state)];
       
       if make_plots > 0
          figure
          plot(theta_grid, ss_cdf);
          title('Theta CDF');
       end
       
       for i = 1:length(quantile_targets)
            if i == 1
                ss_cdf_lb = 0;
                ql = 0;
                agg_prob(i,1) = quantile_targets(i);
            else
                ss_cdf_lb = max(ss_cdf(ss_cdf <= quantile_targets(i-1)));
                % need to handle possibility that first point has a ton of
                % mass on it...
                if isempty(ss_cdf_lb)
                    ss_cdf_lb = 0;
                end
                ql = quantile_targets(i-1);
                agg_prob(i,1) = quantile_targets(i)-ql;
            end    
            ss_cdf_ub = min(ss_cdf(ss_cdf >= quantile_targets(i)));
            if length(ss_cdf_ub) == 0
               ss_cdf_ub = length(ss_cdf); 
            end
            bin_indices = (sum(ss_cdf <= ss_cdf_lb)+1):sum(ss_cdf <= ss_cdf_ub);
           
            % next assign weights associated with each of the intervals.
            % This will be robust to the presence of mass points in the
            % distributions
            if length(bin_indices) == 1
                pweights = 1;
            else
                pweights = steady_state(bin_indices);
                pweights(1) = ss_cdf(bin_indices(1))- ql;
                pweights(end) = quantile_targets(i)-ss_cdf(bin_indices(end-1));
                pweights = pweights / sum(pweights);
            end
        
           % compute average wage growth 
           agg_shock_exposed_wg(i,1) = exposed_shockwg(bin_indices, 1)' * pweights;
           agg_noshock_exposed_wg(i,1) = exposed_cfwg(bin_indices, 1)' * pweights;      
           agg_shock_unexposed_wg(i,1) = unexposed_shockwg(bin_indices, 1)' * pweights;
           agg_noshock_unexposed_wg(i,1) = unexposed_cfwg(bin_indices, 1)' * pweights; 

           % average absolute wage growth
           agg_shock_exposed_awg(i,1) = exposed_shockwg(bin_indices, 2)' * pweights;
           agg_noshock_exposed_awg(i,1) = exposed_cfwg(bin_indices, 2)' * pweights;      
           agg_shock_unexposed_awg(i,1) = unexposed_shockwg(bin_indices, 2)' * pweights;
           agg_noshock_unexposed_awg(i,1) = unexposed_cfwg(bin_indices, 2)' * pweights; 

           % probability of a wage decline larger than p10
           agg_shock_exposed_pctile(i,1) = exposed_shock_lt_pctile_probs(bin_indices, 1)' * pweights;
           agg_noshock_exposed_pctile(i,1) = exposed_cf_lt_pctile_probs(bin_indices, 1)' * pweights;      
           agg_shock_unexposed_pctile(i,1) = unexposed_shock_lt_pctile_probs(bin_indices, 1)' * pweights;
           agg_noshock_unexposed_pctile(i,1) = unexposed_cf_lt_pctile_probs(bin_indices, 1)' * pweights; 

           % expected wage growth & absolute wage growth by income bin
           ss_wg_by_income(i, 1) = sswg(bin_indices, 1)' * pweights;
           ss_awg_by_income(i, 1) = sswg(bin_indices, 2)' * pweights;
       end
      
       
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
       
       diff_coefs = [agg_shock_unexposed_wg - agg_noshock_unexposed_wg, agg_shock_exposed_wg - agg_noshock_exposed_wg].* irf_scale_factor;
       
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
              
       % Compute the E[x * x'] moments needed for the ols formula
       first_col_off_diag_elems = [omega * agg_prob', omega * alpha / p_z * agg_prob']';
       
       xtx = [[omega; first_col_off_diag_elems], [first_col_off_diag_elems'; ...
           [[diag(agg_prob'); ...
           diag(omega * alpha / p_z * agg_prob')],...
           [diag(omega * alpha / p_z * agg_prob'); ...
           diag(omega * alpha / p_z * agg_prob')]]]];
       
       
       wage_growth = (xtx) \ e_xy_wg;
       wage_growth = wage_growth((end - 4):end) .* irf_scale_factor;
       reg_coefs = wage_growth;

       abs_wage_growth = (xtx) \ e_xy_awg;
       abs_wage_growth = abs_wage_growth((end - 4):end) .* irf_scale_factor;
       
       tenth_pctile_probs = (xtx) \ e_xy_pctile;
       tenth_pctile_probs = tenth_pctile_probs((end - 4):end) .* irf_scale_factor;
       
   if make_plots > 0
       
      figure
      
      
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
      
      figure
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
   
   lshare = (H_star * high_wage_ss + L_star * low_wage_ss) / Y_star;
   lshare_shock = (H_shock * high_wage_shock + L_shock * low_wage_shock) / Y_shock;
      
   lshare_irf_sign = ((log(lshare_shock) - log(lshare)) > 0) * abs(log(lshare_shock) - log(lshare)) * 100;
   
   lshare_irf = (log(lshare_shock) - log(lshare)) * agg_scale_factor;
   
   premium = top_twentyfive_wages / bottom_twentyfive_wages;
   
   if make_plots > 0
       
      figure
      theta_plot_x_axis = (high_wage_ss * theta_grid + (low_wage_ss) * (1 - theta_grid) ) ./ ...
          ((low_wage_ss * (1 - min(theta_grid))) + high_wage_ss * min(theta_grid));
      plot(cumsum(steady_state), theta_plot_x_axis, "-o");
      title("Density by Wage / W_l")
      xlabel('Cumulative Density on Theta Grid')
      ylabel('Wages normalized by Low Wage')
      
   end
   
  tenth_pctile_gradient = tenth_pctile_probs(5) - tenth_pctile_probs(1); % changed this from ratio to difference
  wage_gradient = wage_growth(5) - wage_growth(4); % top income bin - second highest income bin gradient target
   
   mom = [lshare; premium; y_irf; lshare_irf; y_irf_sign; lw_irf_sign; ...
       lshare_irf_sign; ...
       abs_wage_growth; ...
       wage_growth; ...
       wage_gradient; ...
       ss_wg_by_income; ...
       ss_awg_by_income; ...
       expected_wage_growth; ...
       expected_abs_wage_growth; ...
       tenth_pctile_probs;
       tenth_pctile_gradient;
       aggregate_sd];
end