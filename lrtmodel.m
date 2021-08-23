function loss = lrtmodel(paramvec, H_inside, make_plots, n_gridpoints,parse_fcn_name)

% transpose for particleswarm
if nargin < 5
    paramvec = paramvec';

    % 0 is geometric increasing from theta
    % 1 is geometric decreasing from 1
    theta_order = 0;

    % set up variables in accordance w/ our model

    % param vec = [phi; alpha; omega; rho; sigma;]

    % rate of moving up the ladder
    phi = (paramvec(1,:));
    % phi = paramvec(1,:);

    % conditional fall probability
    alpha = (paramvec(2,:));
    % alpha = paramvec(2,:);

    lambda = paramvec(9,:);
    mu = paramvec(8,:);

    hc_loss = paramvec(6,:);

    % unconditional growth

    % number of iterations associated with a single shock
    % since a shock is 5 years, this corresponds to there being
    % 4 iterations of the VAR a year
    n_periods = 60;

    g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;

    % death rate
    % corresponds to an average 40 year working life
    delta =  exp(log(0.025 + 0.02 + 1) / (n_periods / 5)) - 1;

    % arrival rate of technology
    % in this case this is also the probability
    % at period t
    % of state 1 happening in period t + 1
    % these shocks are IID so this is true for both
    % initial states
    omegatalpha = (paramvec(3,:));
    omega = omegatalpha / alpha;
    % omega = paramvec(3,:);

    % sigma is outer nest exponent
    % rho is inner nest exponent
    if H_inside == 1
        sigma = (paramvec(5,:));
        rho = sigma - (paramvec(4,:));
    else
        rho = (paramvec(5,:));
        sigma = rho - (paramvec(4,:));
    end

    v = 1;

    % rho = 0.25;
    % sigma = paramvec(5,:);
    % rho = paramvec(4,:);

    % size of technology shock
    % this implies an asymptotic
    % mean of 1 for the xi process

    % TODO: BREAK THIS LINK !!
    kappa = g / omega;


    % kappa = 0.01;
    % kappa = 0.01;

    % number of theta states
    % used to approximate continuous density

    % reversed exponential growth per Dimitris
    % top_grid = - normcdf(paramvec(7,:)) * 5;
    theta0 = 0.05;

    p_z = exp(paramvec(7,:) + log(alpha)) / (1 + exp(paramvec(7,:) + log(alpha)));
    if theta_order == 0
        % growth_rate = paramvec(8,:);
        % n_gridpoints = floor(-log(theta0) / log(1 + growth_rate));
        growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;
        theta_grid = (theta0).*((1 + growth_rate).^(1:(n_gridpoints))); 
    else 
        growth_rate = exp((log(theta0)) / n_gridpoints) - 1;
        theta_grid = 1 - (1.*((1 + growth_rate).^(1:(n_gridpoints))));
    end
    xi_constant = 0;
else
    eval(['[phi,alpha,lambda,mu,hc_loss,n_periods,g,delta,omega,sigma,rho,v,p_z,kappa,theta_grid,theta0,xi_constant] = ', ...
            parse_fcn_name,'(paramvec,H_inside,n_gridpoints);']);
    
end


% need that single obs for xi
n_coefs = 1 + n_gridpoints;

% VAR Intercept Term
% first term corresponds to xi
A_0 = zeros(n_coefs, n_coefs);

% xi depreciation
A_0(1,1) = 1 - g;

% fill in theta section
for i = 2:n_coefs 
   A_0(i, i) = (1 - phi);
   if i < n_coefs 
       A_0(i, i + 1) = phi;
   else
       % can't go past 1
       A_0(i, i) = 1;
   end
end

A_0(2:end, 2:end) = A_0(2:end, 2:end) * (1 - delta);
A_0(2:end,2) = A_0(2:end,2) + delta;



% VAR Intercept Term
% first term corresponds to xi
A_1 = zeros(n_coefs, n_coefs);

% xi depreciation
A_1(1,1) = 1 - g;


% Displacement shock happens last. Below, we will redistribute mass from
% alpha * A_0 across different columns, incorporating the hc loss
A_1(2:end, 2:end) = A_0(2:end, 2:end) * (1 - alpha);



% VAR Intercept Term
% first term corresponds to xi
A_0_no_delta = zeros(n_coefs, n_coefs);

% xi depreciation
A_0_no_delta(1,1) = 1 - g;

% fill in theta section
for i = 2:n_coefs 
   A_0_no_delta(i, i) = (1 - phi);
   if i < n_coefs 
       A_0_no_delta(i, i + 1) = phi;
   else
       % can't go past 1
       A_0_no_delta(i, i) = 1;
   end
end


% VAR Intercept Term
% first term corresponds to xi
A_1_no_delta = zeros(n_coefs, n_coefs);

% xi depreciation
A_1_no_delta(1,1) = 1 - g;

A_1_no_delta(2:end, 2:end) = A_0_no_delta(2:end, 2:end) * (1 - alpha);


% VAR Intercept Term
% first term corresponds to xi
A_0_no_delta_pz = zeros(n_coefs, n_coefs);

% xi depreciation
A_0_no_delta_pz(1,1) = 1 - g;

% fill in theta section
for i = 2:n_coefs 
   A_0_no_delta_pz(i, i) = (1 - phi);
   if i < n_coefs 
       A_0_no_delta_pz(i, i + 1) = phi;
   else
       % can't go past 1
       A_0_no_delta_pz(i, i) = 1;
   end
end


% VAR Intercept Term, for "exposed" workers
% first term corresponds to xi
A_1_no_delta_pz = zeros(n_coefs, n_coefs);

% xi depreciation
A_1_no_delta_pz(1,1) = 1 - g;

A_1_no_delta_pz(2:end, 2:end) = A_0_no_delta(2:end, 2:end) * (1 - p_z);

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
   
   % TODO: this seems to only address the diagonal. What about one above
   % the diagonal (since people learn too)?
   A_1(2:end,upper_fall_index(i,:) + 1) = A_1(2:end,upper_fall_index(i,:) + 1)...
       + alpha * upper_fall_weight(i,:)*A_0(2:end,i+1);
   A_1(2:end,lower_fall_index(i,:) + 1) = A_1(2:end,lower_fall_index(i,:) + 1)...
       + alpha * lower_fall_weight(i,:)*A_0(2:end,i+1);
   
   A_1_no_delta(2:end,upper_fall_index(i,:) + 1) = A_1_no_delta(2:end,upper_fall_index(i,:) + 1)...
       + alpha * upper_fall_weight(i,:)*A_0_no_delta(2:end,i+1);
   A_1_no_delta(2:end,lower_fall_index(i,:) + 1) = A_1_no_delta(2:end,lower_fall_index(i,:) + 1)...
       + alpha * lower_fall_weight(i,:)*A_0_no_delta(2:end,i+1);
   
  A_1_no_delta_pz(2:end,upper_fall_index(i,:) + 1) = A_1_no_delta_pz(2:end,upper_fall_index(i,:) + 1)...
       + p_z * upper_fall_weight(i,:)*A_0_no_delta_pz(2:end,i+1);
   A_1_no_delta_pz(2:end,lower_fall_index(i,:) + 1) = A_1_no_delta_pz(2:end,lower_fall_index(i,:) + 1)...
       + p_z * lower_fall_weight(i,:)*A_0_no_delta_pz(2:end,i+1);
end

% transpose for use w/ Bianchi formulas
% VAR format

A_0 = A_0';
A_1 = A_1';

% transpose for use w/ Bianchi formulas
% VAR format

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

A_0_tilde(2:end,2:end) = A_0_tilde(2:end,2:end)- repmat(c_0_tilde(2:end,1),1,size(c_0_tilde,1)-1);
A_1_tilde(2:end,2:end) = A_1_tilde(2:end,2:end)- repmat(c_1_tilde(2:end,1),1,size(c_1_tilde,1)-1);

% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde_no_delta = A_0_no_delta(1:(end-1),1:(end-1));
A_1_tilde_no_delta = A_1_no_delta(1:(end-1),1:(end-1));

c_0_tilde_no_delta = [xi_constant; A_0_no_delta(2:end-1,end)];
c_1_tilde_no_delta = [kappa+xi_constant; A_1_no_delta(2:end-1,end)];

A_0_tilde_no_delta(2:end,2:end) = A_0_tilde_no_delta(2:end,2:end)- repmat(c_0_tilde_no_delta(2:end,1),1,...
    size(c_0_tilde_no_delta,1)-1);
A_1_tilde_no_delta(2:end,2:end) = A_1_tilde_no_delta(2:end,2:end)- repmat(c_1_tilde_no_delta(2:end,1),1,...
    size(c_1_tilde_no_delta,1)-1);

% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde_no_delta_pz = A_0_no_delta_pz(1:(end-1),1:(end-1));
A_1_tilde_no_delta_pz = A_1_no_delta_pz(1:(end-1),1:(end-1));

c_0_tilde_no_delta_pz = [xi_constant; A_0_no_delta_pz(2:end-1,end)];
c_1_tilde_no_delta_pz = [kappa+xi_constant; A_1_no_delta_pz(2:end-1,end)];

A_0_tilde_no_delta_pz(2:end,2:end) = A_0_tilde_no_delta_pz(2:end,2:end)- ...
    repmat(c_0_tilde_no_delta_pz(2:end,1),1,size(c_0_tilde_no_delta_pz,1)-1);
A_1_tilde_no_delta_pz(2:end,2:end) = A_1_tilde_no_delta_pz(2:end,2:end)- ...
    repmat(c_1_tilde_no_delta_pz(2:end,1),1,size(c_1_tilde_no_delta_pz,1)-1);

% transition matrix
T = [[1 - omega, omega]; [1 - omega, omega]];

% stationary distribution
[V,D] = eigs(T');
DVec  = diag(D);
piVec = V(:, abs(DVec-1) == min(abs(DVec-1)));
piVec = piVec ./ sum(piVec);

H = T';

C = blkdiag(c_0_tilde, c_1_tilde);
%C = blkdiag([0, 0, repelem(0, n_gridpoints - 1)]', [kappa, 0, repelem(0, n_gridpoints - 1)]');

bianchi_omega = blkdiag(A_0_tilde, A_1_tilde) * kron(H,eye(n_coefs-1));
q = (eye((n_coefs-1)*2) - bianchi_omega) \ (C * piVec);         % eq. (3)
bianchi_omegatilde = [bianchi_omega, C*H; zeros(2, (n_coefs-1)*2), H];  % eq. (5)

w = repmat(eye((n_coefs-1)),1,2);
mu_ss = w * q;

wtilde = [w, zeros((n_coefs-1),2)];
qtilde = [q; piVec];

steady_state = [mu_ss(2:(n_coefs-1)); 1 - sum(mu_ss(2:end))];
figure(1)
plot(theta_grid,steady_state)
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
%     kappa, rho, sigma, alpha, phi, xi_var, ...
%     A_0_tilde, A_1_tilde, ...
%     c_0_tilde, c_1_tilde, omega, n_periods, v, ...
%     A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
%     c_1_tilde_no_delta, A_0_tilde_no_delta_pz, ...
%     A_1_tilde_no_delta_pz, c_0_tilde_no_delta_pz, ...
%     c_1_tilde_no_delta_pz, p_z), [0; 0], options); 
% normcdf(lambdamu)

theor_mom = calcmom(lambda, mu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, ...
    A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta, A_0_tilde_no_delta_pz, ...
    A_1_tilde_no_delta_pz, c_0_tilde_no_delta_pz, ...
    c_1_tilde_no_delta_pz, p_z, 1, make_plots, A_1);



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
bottom_density_loss = (steady_state(1) > 0.1) * abs((steady_state(1) - 0.1)) * 100;
% half_income_from_low_skill = low_wage *(1-theta_0) / (low_wage * (1-theta_0) +  theta_0 * high_wage) >= 0.5;


expected_wage_growth_by_income = [0.006216; -0.08353; -0.08933; -0.09111; -0.1197];
expected_abs_wage_growth_by_income = [0.4466; 0.2961; 0.2559; 0.2437; 0.2792];

emp_mom = [0.66; 2.45; 0.0281; -0.0125; 0; 0; 0; ...
             emp_abs_wage_growth; ...
             emp_wage_growth; ...
             expected_wage_growth_by_income; ...
             expected_abs_wage_growth_by_income; ...
             -0.06313; 0.3171; ... expected wage growth, expected abs wage growth
             tenth_pctile_probs];
         
weight_vec = [20; 10; 25; 25; 1; 1; 1;
         5; 5; 5; 5; 5; ...
         15; 10; 10; 10; 15; ...
         5; 5; 5; 5; 5; ...
         5; 5; 5; 5; 5; ...
         10; 10; ...
         0; 0; 0; 0; 0];

loss_vec = (theor_mom - emp_mom) ./ (0.01 + abs(emp_mom)) .* weight_vec;
%  bars(labels, loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
%  
% [emp_mom, theor_mom]

if make_plots > 0
   figure(2)
   momlabels = categorical(1:31, 1:31, {'Labor Share', 'Wage Ratios', 'Output IRF','LShare IRF',...
          'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
           'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]',...
           'E(WG[0,25])','E(WG[25,50])','E(WG[50,75])','E(WG[75,95])','E(WG[95,100])',...
           'E(AWG[0,25])','E(AWG[25,50])','E(AWG[50,75])','E(AWG[75,95])','E(AWG[95,100])', ...
           'E(WG)', 'E(AWG)', ...
           'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]'},...
           'Ordinal',true);
      bar(momlabels(3:end)', [theor_mom([3:4, 8:(end)]), emp_mom([3:4, 8:(end)])])
     title('Moment Matching (excluding signs & levels)')
     
     figure(1)
     bar(momlabels(1:2)', [theor_mom(1:2), emp_mom(1:2)])
     title('Labor Share & Wage Ratio')
     
     figure(3)
     bar(momlabels(1:end)', weight_vec([1:4, 8:(end)]).^2)
     title('Weights)')
     
     
     
   % scaling factors to convert from one shock units to 1 SD units
   agg_scale_factor = sqrt(n_periods / 5) * sqrt(omega * (1 - omega));
   irf_scale_factor = sqrt(n_periods / 5) * sqrt(omega * alpha / p_z * (1 - omega * alpha / p_z));
     
    names = {'HC Increase Prob', 'Conditional Fall Prob', ...
        'Shock Prob', 'Skilled Share', 'Technology Share', 'Skilled Curvature', ...
        'Unskilled Curvature', 'DRS Param', 'Bottom Rung', 'P(fall | shock, exposed)',...
        'kappa','Xi intercept','Xi shock size (annualized)','Xi mean','Xi std dev','Human capital loss'};

    all_params = [phi, alpha, omega, mu, lambda, sigma, rho, v, theta0, p_z, ...
                  kappa,xi_constant,kappa*agg_scale_factor, xi_star, sqrt(xi_var), hc_loss]';

    disp(table(names', all_params))
end 
     
loss_vec = [loss_vec; bottom_density_loss; top_density_loss];

% miss = ([theor_mom([1:2, 6:end]) - emp_mom([1:2, 6:end])] ./ (0.01 + abs(emp_mom([1:2, 6:end])))).^2;
% miss = miss .* weight_vec ./ sum(miss .* weight_vec);
% 
% figure(3)
if make_plots > 0
labels = categorical(1:36, 1:36, {'Labor Share', 'Wage Ratio','Output IRF','LShare IRF',...
    'Output ','Wage Sign', 'Lshare IRF sign', ...
     'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
     'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]',...
     'E(WG[0,25])','E(WG[25,50])','E(WG[50,75])','E(WG[75,95])','E(WG[95,100])',...
     'E(AWG[0,25])','E(AWG[25,50])','E(AWG[50,75])','E(AWG[75,95])','E(AWG[95,100])', ...
     'E(WG)', 'E(AWG)', ...
     'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]',...
     'Bottom Density Penalty', 'Top Density Penalty'}, 'Ordinal',true);

figure(4)
bar(labels', loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
title('Weighted Percent Loss Contribution')
end
if any(isnan(loss_vec) | ~isreal(loss_vec))
    loss = 1e16;
else
    loss = loss_vec' * loss_vec;
end

