function [loss, emp_mom, theor_mom, wh, wl, Y, H, lshare, xi, twenty_fifth_pctile, seventy_fifth_pctile, ...
percent_from_top_five, percent_from_top_one, ss_dist, ss_wage_grid] =...
    lrtmodel_shock_omega_phi(paramvec, make_plots, hyperparams, omega_shock_amount, phi_shock_amount)

H_inside = hyperparams.H_inside;
n_gridpoints = hyperparams.n_gridpoints; 
scale_period = hyperparams.scale_period;
n_periods = hyperparams.n_periods;
parse_fcn_name = hyperparams.parse_fcn_name;

if isempty(parse_fcn_name)
   [phi,alpha_param,lambda,mu,hc_loss,n_periods,g,delta,omega,sigma_param,rho,v,p_z,kappa,theta_grid,theta0,xi_constant, p0_share, p_up, p_down] = ...
        parse_model_params_v4(paramvec, hyperparams);
else    
eval(['[phi,alpha_param,lambda,mu,hc_loss,n_periods,g,delta,omega,sigma_param,rho,v,p_z,kappa,theta_grid,theta0,xi_constant, p0_share, p_up, p_down] = ', ...
        parse_fcn_name,'(paramvec, hyperparams);']);
end 

% need that single obs for xi
n_coefs = 2 + n_gridpoints;

% VAR Intercept Term
% first term corresponds to xi
A_0 = zeros(n_coefs, n_coefs);

% xi depreciation
A_0(1,1) = 1 - g;
A_0(2, 2) = 1;

% fill in theta section
% phi is probability of moving
% p_up is probability that the move is upwards
% p_down is probability that the move is downwards
for i = 3:n_coefs 
   A_0(i, i) = (1 - phi);
   if i == 3
     A_0(i, i) = 1 - phi * p_up;
     A_0(i, i + 1) = phi * p_up;

   elseif i < n_coefs 
     A_0(i, i + 1) = phi * p_up;
     A_0(i, i - 1) = phi * p_down;
   else
       % can't go past 1
     A_0(i, i - 1) = phi * p_down;
     A_0(i, i) = 1 - phi * p_down;
   end
end

A_0(3:end, 3:end) = A_0(3:end, 3:end) * (1 - delta);
A_0(3:end,3) = A_0(3:end,3) + delta;


% VAR Intercept Term
% first term corresponds to xi
A_1 = zeros(n_coefs, n_coefs);

% xi depreciation
A_1(1,1) = 1 - g;

% folks who don't move
A_1(2,2) = 1;

% Displacement shock happens last. Below, we will redistribute mass from
% alpha_param * A_0 across different columns, incorporating the hc loss
A_1(3:end, 3:end) = A_0(3:end, 3:end) * (1 - alpha_param);

% VAR Intercept Term
% first term corresponds to xi
A_0_no_delta = zeros(n_coefs, n_coefs);

% xi depreciation
A_0_no_delta(1,1) = 1 - g;

% p0 doesn't move
A_0_no_delta(2, 2) = 1;


% fill in theta section
% phi is probability of moving
% p_up is probability that the move is upwards
% p_down is probability that the move is downwards
for i = 3:n_coefs 
   A_0_no_delta(i, i) = (1 - phi);
   if i == 3
     A_0_no_delta(i, i) = 1 - phi * p_up;
     A_0_no_delta(i, i + 1) = phi * p_up;

   elseif i < n_coefs 
     A_0_no_delta(i, i + 1) = phi * p_up;
     A_0_no_delta(i, i - 1) = phi * p_down;
   else
       % can't go past 1
     A_0_no_delta(i, i - 1) = phi * p_down;
     A_0_no_delta(i, i) = 1 - phi * p_down;
   end
end

% VAR Intercept Term
% first term corresponds to xi
A_1_no_delta = zeros(n_coefs, n_coefs);

% xi depreciation
A_1_no_delta(1,1) = 1 - g;

% folks who never move
A_1_no_delta(2,2) = 1;

A_1_no_delta(3:end, 3:end) = A_0_no_delta(3:end, 3:end) * (1 - alpha_param);


% VAR Intercept Term
% first term corresponds to xi
A_0_no_delta_pz = zeros(n_coefs, n_coefs);

% xi depreciation
A_0_no_delta_pz(1,1) = 1 - g;

% folks who don't move
A_0_no_delta_pz(2,2) = 1;

% fill in theta section
% phi is probability of moving
% p_up is probability that the move is upwards
% p_down is probability that the move is downwards
for i = 3:n_coefs 
   A_0_no_delta_pz(i, i) = (1 - phi);
   if i == 3
     A_0_no_delta_pz(i, i) = 1 - phi * p_up;
     A_0_no_delta_pz(i, i + 1) = phi * p_up;

   elseif i < n_coefs 
     A_0_no_delta_pz(i, i + 1) = phi * p_up;
     A_0_no_delta_pz(i, i - 1) = phi * p_down;
   else
       % can't go past 1
     A_0_no_delta_pz(i, i - 1) = phi * p_down;
     A_0_no_delta_pz(i, i) = 1 - phi * p_down;
   end
end


% VAR Intercept Term, for "exposed" workers
% first term corresponds to xi
A_1_no_delta_pz = zeros(n_coefs, n_coefs);

% xi depreciation
A_1_no_delta_pz(1,1) = 1 - g;
 
% folks who don't move
A_1_no_delta_pz(2, 2) = 1;

A_1_no_delta_pz(3:end, 3:end) = A_0_no_delta_pz(3:end, 3:end) * (1 - p_z);

% transpose for use w/ Bianchi formulas
% VAR format

% TODO: probably better to make this exp(-hc_loss) instead...
new_theta = theta_grid - theta_grid * hc_loss;

lower_fall_index = zeros(size(theta_grid))';
upper_fall_index = zeros(size(theta_grid))';
lower_fall_weight = zeros(size(theta_grid))';
upper_fall_weight = zeros(size(theta_grid))';

for i = 1:length(new_theta) 
   new_gridpoint_l = max(cumsum(theta_grid < new_theta(i)));
   if new_gridpoint_l == 0
      lower_fall_index(i,:) = 1;
      upper_fall_index(i,:) = 1;
      lower_fall_weight(i,:) = 1;
      upper_fall_weight(i,:) = 0;
   else 
      new_gridpoint_u = new_gridpoint_l + 1;
      lower_fall_index(i,:) = new_gridpoint_l;
      upper_fall_index(i,:) = new_gridpoint_u;
      
      upper_fall_weight(i,:) = (log(theta_grid(new_gridpoint_u)) - log(new_theta(i)))...
          / (log(theta_grid(new_gridpoint_u)) - log(theta_grid(new_gridpoint_l)));
      lower_fall_weight(i,:) = 1 - upper_fall_weight(i,:);
   end
   
   A_1(3:end,upper_fall_index(i,:) + 2) = A_1(3:end,upper_fall_index(i,:) + 2)...
       + alpha_param * upper_fall_weight(i,:)*A_0(3:end,i+2);
   A_1(3:end,lower_fall_index(i,:) + 2) = A_1(3:end,lower_fall_index(i,:) + 2)...
       + alpha_param * lower_fall_weight(i,:)*A_0(3:end,i+2);
   
   A_1_no_delta(3:end,upper_fall_index(i,:) + 2) = A_1_no_delta(3:end,upper_fall_index(i,:) + 2)...
       + alpha_param * upper_fall_weight(i,:)*A_0_no_delta(3:end,i+2);
   A_1_no_delta(3:end,lower_fall_index(i,:) + 2) = A_1_no_delta(3:end,lower_fall_index(i,:) + 2)...
       + alpha_param * lower_fall_weight(i,:)*A_0_no_delta(3:end,i+2);
   
  A_1_no_delta_pz(3:end,upper_fall_index(i,:) + 2) = A_1_no_delta_pz(3:end,upper_fall_index(i,:) + 2)...
       + p_z * upper_fall_weight(i,:)*A_0_no_delta_pz(3:end,i+2);
  A_1_no_delta_pz(3:end,lower_fall_index(i,:) + 2) = A_1_no_delta_pz(3:end,lower_fall_index(i,:) + 2)...
       + p_z * lower_fall_weight(i,:)*A_0_no_delta_pz(3:end,i+2);
end

% transpose for use w/ Bianchi formulas
% VAR format

A_0 = A_0';
A_1 = A_1';

A_0_no_delta = A_0_no_delta';
A_1_no_delta = A_1_no_delta';

A_0_no_delta_pz = A_0_no_delta_pz';
A_1_no_delta_pz = A_1_no_delta_pz';


% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde = A_0(1:(end-1),1:(end-1));
A_1_tilde = A_1(1:(end-1),1:(end-1));

c_0_tilde = [xi_constant; A_0(2:end-1,end)];
c_1_tilde = [kappa+xi_constant; A_1(2:end-1,end)];

A_0_tilde(3:end,3:end) = A_0_tilde(3:end,3:end)- repmat(c_0_tilde(3:end,1),1,size(c_0_tilde,1)-2);
A_1_tilde(3:end,3:end) = A_1_tilde(3:end,3:end)- repmat(c_1_tilde(3:end,1),1,size(c_1_tilde,1)-2);

% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde_no_delta = A_0_no_delta(1:(end-1),1:(end-1));
A_1_tilde_no_delta = A_1_no_delta(1:(end-1),1:(end-1));

c_0_tilde_no_delta = [xi_constant; A_0_no_delta(2:end-1,end)];
c_1_tilde_no_delta = [kappa+xi_constant; A_1_no_delta(2:end-1,end)];

A_0_tilde_no_delta(3:end,3:end) = A_0_tilde_no_delta(3:end,3:end)- repmat(c_0_tilde_no_delta(3:end,1),1,...
    size(c_0_tilde_no_delta,1)-2);
A_1_tilde_no_delta(3:end,3:end) = A_1_tilde_no_delta(3:end,3:end)- repmat(c_1_tilde_no_delta(3:end,1),1,...
    size(c_1_tilde_no_delta,1)-2);

% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde_no_delta_pz = A_0_no_delta_pz(1:(end-1),1:(end-1));
A_1_tilde_no_delta_pz = A_1_no_delta_pz(1:(end-1),1:(end-1));

c_0_tilde_no_delta_pz = [xi_constant; A_0_no_delta_pz(2:end-1,end)];
c_1_tilde_no_delta_pz = [kappa+xi_constant; A_1_no_delta_pz(2:end-1,end)];

A_0_tilde_no_delta_pz(3:end,3:end) = A_0_tilde_no_delta_pz(3:end,3:end)- ...
    repmat(c_0_tilde_no_delta_pz(3:end,1),1,size(c_0_tilde_no_delta_pz,1)-2);
A_1_tilde_no_delta_pz(3:end,3:end) = A_1_tilde_no_delta_pz(3:end,3:end)- ...
    repmat(c_1_tilde_no_delta_pz(3:end,1),1,size(c_1_tilde_no_delta_pz,1)-2);

% transition matrix
T = [[1 - omega, omega]; [1 - omega, omega]];

% stationary distribution
[V,D] = eigs(T');
DVec  = diag(D);
piVec = V(:, abs(DVec-1) == min(abs(DVec-1)));
piVec = piVec ./ sum(piVec);

H = T';

% Note: because of the presence of the absorbing state in the first entry
% of the state vector, the solution for the steady state is indeterminate.
% So, we need to make some additional restrictions (drop first elemnt of probability vector)
c_0_tilde2 = c_0_tilde([1, 3:end]);
c_1_tilde2 = c_1_tilde([1, 3:end]);
% need to change the intercepts to scale down by mass in the absorbing
% state, since modeled probabilities need to sum to 1-p0_share, not 1
c_0_tilde2(2:end) = c_0_tilde2(2:end)*(1-p0_share);
c_1_tilde2(2:end) = c_1_tilde2(2:end)*(1-p0_share);

A_0_tilde2 = A_0_tilde([1, 3:end],[1, 3:end]);
A_1_tilde2 = A_1_tilde([1, 3:end],[1, 3:end]);

C = blkdiag(c_0_tilde2, c_1_tilde2);
%C = blkdiag([0, 0, repelem(0, n_gridpoints - 1)]', [kappa, 0, repelem(0, n_gridpoints - 1)]');

bianchi_omega = blkdiag(A_0_tilde2, A_1_tilde2) * kron(H,eye(n_coefs-2));
q = (eye((n_coefs-2)*2) - bianchi_omega) \ (C * piVec);         % eq. (3)
bianchi_omegatilde = [bianchi_omega, C*H; zeros(2, (n_coefs-2)*2), H];  % eq. (5)

w = repmat(eye((n_coefs-2)),1,2);
mu_ss = w * q;

steady_state = [p0_share; [mu_ss(2:end); (1 - p0_share) - sum(mu_ss(2:end))] ];
ss_dist = steady_state;

theta_grid = [0, theta_grid];

if make_plots > 0
    figure
    plot(theta_grid,steady_state)
    title('Steady state theta density')
end
% steady state values
H_star = theta_grid * steady_state;
L_star = 1 - H_star;
xi_star = mu_ss(1);

% figure(2)
% plot(shock_state)
xi_var = kappa^2 / (2 * g - g^2) * (1 - omega) * (omega);

% fsolve() for 2 x 2 system
% choose mu and lambda to match labor share = 0.66
% and ratio of high/low wages
% options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
% lambdamu = fsolve(@(x) calcloss(x, theta_grid, steady_state, xi_star, ...
%     kappa, rho, sigma_param, alpha_param, phi, xi_var, ...
%     A_0_tilde, A_1_tilde, ...
%     c_0_tilde, c_1_tilde, omega, n_periods, v, ...
%     A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
%     c_1_tilde_no_delta, A_0_tilde_no_delta_pz, ...
%     A_1_tilde_no_delta_pz, c_0_tilde_no_delta_pz, ...
%     c_1_tilde_no_delta_pz, p_z), [0; 0], options); 
% normcdf(lambdamu)
theor_mom = calcmom(lambda, mu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma_param, alpha_param, phi, xi_var, ...
    A_0, A_1, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_no_delta, A_1_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta, A_0_no_delta_pz, ...
    A_1_no_delta_pz, c_0_tilde_no_delta_pz, ...
    c_1_tilde_no_delta_pz, p_z, 1, make_plots, scale_period, H_inside);



% we target labor share above, now we target other xistuff
theor_mom = theor_mom(1:end,:);

emp_abs_wage_growth = ... 
    [-0.002551; ...
     0.0009579; ...
     0.004768; ...
     0.007638; ...
     0.02014];

emp_wage_growth = [-0.01486; -0.01008; -0.01178; -0.01167; -0.02467];

tenth_pctile_probs = [0.00286; 0.002619; 0.003889; 0.003941; 0.01255];

top_density_loss = (steady_state(end) > 0.01) * abs((steady_state(end) - 0.01)) * 100;
bottom_density_loss = (steady_state(1) > 0.1) * abs((steady_state(1) - 0.1)) * 00;
% half_income_from_low_skill = low_wage *(1-theta_0) / (low_wage * (1-theta_0) +  theta_0 * high_wage) >= 0.5;


expected_wage_growth_by_income = [0.006216; -0.08353; -0.08933; -0.09111; -0.1197];
expected_abs_wage_growth_by_income = [0.4466; 0.2961; 0.2559; 0.2437; 0.2792];

emp_mom = [0.66; 2.45; 0.0281; -0.0125; 0; 0; 0; ...
             emp_abs_wage_growth; ...
             emp_wage_growth; ...
             -0.0130; ... gradient between 5th and 4th income bin in the data
             expected_wage_growth_by_income; ...
             expected_abs_wage_growth_by_income; ...
             -0.06313; 0.3171; ... expected wage growth, expected abs wage growth
             tenth_pctile_probs; ...
             tenth_pctile_probs(5)-tenth_pctile_probs(1); ... % Difference higher high income p10 vs. lowest income p10
             0.592 / sqrt(60)]; % aggregate standard deviation / sqrt(60)
    
weight_vec = hyperparams.weight_vec;

loss_vec = (theor_mom - emp_mom) ./ (0.01 + abs(emp_mom)) .* weight_vec;
%  bars(labels, loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
%  
% [emp_mom, theor_mom]

omega = omega + omega_shock_amount;
phi = phi + phi_shock_amount;


% VAR Intercept Term
% first term corresponds to xi
A_0_shock_params = zeros(n_coefs, n_coefs);

% xi depreciation
A_0_shock_params(1,1) = 1 - g;
A_0_shock_params(2, 2) = 1;

% fill in theta section
% phi is probability of moving
% p_up is probability that the move is upwards
% p_down is probability that the move is downwards
for i = 3:n_coefs 
   A_0_shock_params(i, i) = (1 - phi);
   if i == 3
     A_0_shock_params(i, i) = 1 - phi * p_up;
     A_0_shock_params(i, i + 1) = phi * p_up;

   elseif i < n_coefs 
     A_0_shock_params(i, i + 1) = phi * p_up;
     A_0_shock_params(i, i - 1) = phi * p_down;
   else
       % can't go past 1
     A_0_shock_params(i, i - 1) = phi * p_down;
     A_0_shock_params(i, i) = 1 - phi * p_down;
   end
end

A_0_shock_params(3:end, 3:end) = A_0_shock_params(3:end, 3:end) * (1 - delta);
A_0_shock_params(3:end,3) = A_0_shock_params(3:end,3) + delta;


% VAR Intercept Term
% first term corresponds to xi
A_1_shock_params = zeros(n_coefs, n_coefs);

% xi depreciation
A_1_shock_params(1,1) = 1 - g;

% folks who don't move
A_1_shock_params(2,2) = 1;


% Displacement shock happens last. Below, we will redistribute mass from
% alpha_param * A_0 across different columns, incorporating the hc loss
A_1_shock_params(3:end, 3:end) = A_0_shock_params(3:end, 3:end) * (1 - alpha_param);

for i = 1:length(new_theta)    
   % TODO: this seems to only address the diagonal. What about one above
   % the diagonal (since people learn too)?
   A_1_shock_params(3:end,upper_fall_index(i,:) + 2) = A_1_shock_params(3:end,upper_fall_index(i,:) + 2)...
       + alpha_param * upper_fall_weight(i,:)*A_0_shock_params(3:end,i+2);
   A_1_shock_params(3:end,lower_fall_index(i,:) + 2) = A_1_shock_params(3:end,lower_fall_index(i,:) + 2)...
       + alpha_param * lower_fall_weight(i,:)*A_0_shock_params(3:end,i+2);

end

% transpose for use w/ Bianchi formulas
% VAR format

A_0_shock_params = A_0_shock_params';
A_1_shock_params = A_1_shock_params';

if make_plots > 0
   figure
   momlabels = categorical(1:34, 1:34, {'Labor Share', 'Wage Ratios', 'Output IRF','LShare IRF',...
          'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
           'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]', 'WG[95,100] - WG[75,95]', ...
           'E(WG[0,25])','E(WG[25,50])','E(WG[50,75])','E(WG[75,95])','E(WG[95,100])',...
           'E(AWG[0,25])','E(AWG[25,50])','E(AWG[50,75])','E(AWG[75,95])','E(AWG[95,100])', ...
           'E(WG)', 'E(AWG)', ...
           'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]', 'P10 Gradient', 'Aggregate SD'},...
           'Ordinal',true);
      bar(momlabels([3:15, (end - 6):(end)])', [theor_mom([3:4, 8:(18),(end - 6):(end)]), ...
          emp_mom([3:4, 8:(18), (end - 6):(end)])])
     title('Moment Matching (excluding signs & levels)')
     
     figure
     bar(momlabels([1:2])', [theor_mom([1:2]), emp_mom([1:2])])
     title('Labor Share & Wage Ratio')
     
     figure
     bar(momlabels(1:end)', weight_vec([1:4, 8:end]).^2)
     title('Weights)')
     
     
     
   % scaling factors to convert from one shock units to 1 SD units
   agg_scale_factor = sqrt(scale_period) * sqrt(omega * (1 - omega));
   irf_scale_factor = sqrt(scale_period) * sqrt(omega * alpha_param / p_z * (1 - omega * alpha_param / p_z));
     
    names = {'HC Increase Prob', 'Conditional Fall Prob', ...
        'Shock Prob', 'Skilled Share', 'Technology Share', 'Skilled Curvature', ...
        'Unskilled Curvature', 'DRS Param', 'Bottom Rung', 'P(fall | shock, exposed)',...
        'kappa','Xi intercept','Xi shock size (annualized)','Xi mean','Xi std dev',...
        'Human capital loss', 'g', 'Bottom Rung Share','P_up'};

    all_params = [phi, alpha_param, omega, mu, lambda, sigma_param, rho, v, theta0, p_z, ...
                  kappa,xi_constant,kappa*agg_scale_factor, xi_star, sqrt(xi_var), hc_loss, g, p0_share,p_up]';

    disp(table(names', all_params))
end 
     
loss_vec = [loss_vec; bottom_density_loss; top_density_loss];

% miss = ([theor_mom([1:2, 6:end]) - emp_mom([1:2, 6:end])] ./ (0.01 + abs(emp_mom([1:2, 6:end])))).^2;
% miss = miss .* weight_vec ./ sum(miss .* weight_vec);
% 
% figure(3)


   dist_over_time = steady_state;

   H_star = theta_grid * steady_state;
   L_star = 1 - H_star;
   
   % TODO: should we incorporate the constants for kappa more explicitly? 
   shock_vars = A_1 * [xi_star; steady_state];
   
   xi_shock = shock_vars(1) + c_1_tilde(1);
   shock_state = shock_vars(2:end);

   % steady state production values
   X_star = calc_X(xi_star, H_star, L_star, lambda, rho, H_inside);
   Y_star = calc_Y(H_star, L_star, X_star, mu, sigma_param, v, H_inside);

   % high and low wages at steady state
   high_wage = w_h(H_star, L_star, xi_star, rho, sigma_param, mu, lambda, v, H_inside);
   low_wage = w_l(H_star, L_star, xi_star, rho, sigma_param, mu, lambda, v, H_inside);
   H = zeros(n_periods + 1, 1);
   L = zeros(n_periods + 1, 1);
   xi = zeros(n_periods + 1, 1);
   
   wh = zeros(n_periods + 1, 1);
   wl = zeros(n_periods + 1, 1);
   
   X = zeros(n_periods + 1, 1);
   Y = zeros(n_periods + 1, 1);
   
   
   lshare = zeros(n_periods + 1, 1);
   
   H(1) = H_star;
   L(1) = L_star;
   xi(1) = xi_star;
   wh(1) = high_wage;
   wl(1) = low_wage;
   Y(1) = Y_star;
   X(1) = X_star;
   lshare(1) = (H(1) * wh(1) + L(1) * wl(1)) / Y_star;
   
   shock_vec = [xi_star; steady_state(1:(end))];
   
   % this just adds kappa in times of trouble
   % and the xi intercept term in good times
   int_for_xi_0 = zeros(size(A_0, 1), 1);
   int_for_xi_0(1) = c_0_tilde(1);
   int_for_xi_1 = zeros(size(A_0, 1), 1);
   int_for_xi_1(1) = c_1_tilde(1);
   
   
   for i=1:(scale_period * 60)
       shock_vec = (1 - omega) * (A_0_shock_params * shock_vec + int_for_xi_0) + ...
           omega  * (A_1_shock_params * shock_vec + int_for_xi_1);
          
       xi(i + 1) = shock_vec(1,:);
       shock_state = shock_vec(3:end,:);
       dist_over_time = [dist_over_time, shock_vec(2:end)];

       shock_vec_for_calcs = [p0_share; shock_state];

       H(i + 1) = theta_grid * shock_vec_for_calcs;
       L(i + 1) = 1 - H(i + 1 );

       % shock state production values
       X(i + 1) = calc_X(xi(i + 1),H(i + 1), L(i + 1), lambda, rho, H_inside);
       Y(i + 1) = calc_Y(H(i + 1), L(i + 1), X(i + 1), mu, sigma_param, v, H_inside);


       % high and low wages at shock state
       wh(i + 1) = w_h(H(i + 1), L(i + 1), xi(i + 1), rho, sigma_param, mu, lambda, v, H_inside);
       wl(i + 1) = w_l(H(i + 1), L(i + 1), xi(i + 1), rho, sigma_param, mu, lambda, v, H_inside);
       lshare(i + 1) = (H(i + 1) * wh(i + 1) + L(i + 1) * wl(i + 1)) / Y(i + 1);
   end
   
    wages = wh * theta_grid + wl * (1 - theta_grid);

    wage_ts = zeros(size(dist_over_time, 2), 1);
    top_five_wages = zeros(size(dist_over_time, 2), 1);
    top_one_wages = zeros(size(dist_over_time, 2), 1);
    for t = 1:size(dist_over_time, 2)
       wage_ts(t, :) = wages(t, :) * dist_over_time(:, t);

       
       fifth_pts = cumsum(dist_over_time(:, t)) > 0.95;
       first_pts = cumsum(dist_over_time(:, t)) > 0.99;
       
       if sum(dist_over_time(fifth_pts, t)) < 0.05
           fifth_pts(sum(fifth_pts == 0)) = 1;
       end
       
       if sum(dist_over_time(first_pts, t)) < 0.01
            first_pts(sum(first_pts == 0)) = 1;
       end
       
       
       fifth_dist = dist_over_time(fifth_pts, t);
       fifth_dist(1) = fifth_dist(1) - (sum(fifth_dist) - 0.05);
       
       first_dist = dist_over_time(first_pts, t);
       first_dist(1) = first_dist(1) - (sum(first_dist) - 0.01);
       
       top_five_wages(t, :) = wages(t, fifth_pts) * fifth_dist;
       top_one_wages(t, :) = wages(t,first_pts) * first_dist;

    end

    percent_from_top_one = top_one_wages ./ wage_ts;
    percent_from_top_five = top_five_wages ./ wage_ts;
   
    twenty_fifth_pctile = zeros(size(dist_over_time, 2), 1);
    seventy_fifth_pctile = zeros(size(dist_over_time, 2), 1);
    for i=1:(scale_period * 60)
        current_wages = wh(i) * theta_grid + wl(i) * (1 - theta_grid);
        twenty_fifth_pctile(i,:) = wl(i);
        seventy_fifth_pctile(i,:) =  weighted_quantile_between_grid(current_wages, dist_over_time(:, i), 0.75);
    end
    
   if make_plots > 0
    labels = categorical(1:39, 1:39, {'Labor Share', 'Wage Ratio','Output IRF','LShare IRF',...
        'Output ','Wage Sign', 'Lshare IRF sign', ...
         'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
         'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]', 'WG[95,100] - WG[75,95]',...
         'E(WG[0,25])','E(WG[25,50])','E(WG[50,75])','E(WG[75,95])','E(WG[95,100])',...
         'E(AWG[0,25])','E(AWG[25,50])','E(AWG[50,75])','E(AWG[75,95])','E(AWG[95,100])', ...
         'E(WG)', 'E(AWG)', ...
         'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]', 'P10 Gradient', 'Aggregate SD',...
         'Bottom Density Penalty', 'Top Density Penalty'}, 'Ordinal',true);

    figure
    bar(labels', loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
    title('Weighted Percent Loss Contribution')
    addpath('matlab2tikz/src/')

    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (wh(1:scale_period * 60)./wh(1) - 1)*agg_scale_factor, '.-')
    title("High Wage")

    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (wl(1:scale_period * 60)./wl(1) - 1)*agg_scale_factor, '.-')
    title("Low Wage")


    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (L(1:scale_period * 60)./L(1) - 1)*agg_scale_factor, '.-')
    title("L Skill Level")

    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (H(1:scale_period * 60)./H(1) - 1)*agg_scale_factor, '.-')
    title("H Skill Level")

    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (xi(1:scale_period * 60)./xi(1) - 1)*agg_scale_factor, '.-')
    title("Technology Level")

    wage_diff = wh - wl;
    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (wage_diff(1:scale_period * 60)./wage_diff(1) - 1)*agg_scale_factor, '.-')
    title("High Wage - Low Wage")

    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (lshare(1:scale_period * 60)./lshare(1) - 1)*agg_scale_factor, '.-')
    title('Labor Share')


    figure
    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (X(1:scale_period * 60)./X(1) - 1)*agg_scale_factor, '.-')
    title("Composite Good")

    plot([0, (1:(scale_period * 60 - 1))]'./scale_period, (Y(1:scale_period * 60)./Y(1) - 1)*agg_scale_factor, '.-')
    title("Output Level")

    plot([0, (1:(scale_period * 60))]'./scale_period, [twenty_fifth_pctile, seventy_fifth_pctile], '.-')
    title("Inequality over Time")

    level_names = {'High Wage', 'Low Wage', 'Wage Difference', ...
        'Output', 'Composite Good', 'Low Skill Amount', 'High Skill', ...
        'Labor Share',  'Xi', 'Share of Wages from L'};

    all_levels = [wh(1), wl(1), wage_diff(1), Y(1), X(1), L(1), H(1), lshare(1), xi(1), wl(1) * L(1) / (wl(1) * L(1) + wh(1) * H(1))]';

    disp(table(level_names', all_levels))



   end

if any(isnan(loss_vec) | ~isreal(loss_vec))
    loss = 1e16;
else
    loss = loss_vec' * loss_vec;
end

