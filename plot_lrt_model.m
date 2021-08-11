
% transpose for particleswarm
paramvec = paramvec';

% set up variables in accordance w/ our model

% param vec = [phi; alpha; omega; rho; sigma;]

% rate of moving up the ladder
phi = (paramvec(1,:));
% phi = paramvec(1,:);

% conditional fall probability
alpha = (paramvec(2,:));
% alpha = paramvec(2,:);

% unconditional growth

% number of iterations associated with a single shock
% since a shock is 5 years, this corresponds to there being
% 4 iterations of the VAR a year
n_periods = 20;

g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;

% death rate
% corresponds to an average 40 year working life
delta =  exp(log(0.025 + 1) / (n_periods / 5)) - 1;

% arrival rate of technology
% in this case this is also the probability
% at period t
% of state 1 happening in period t + 1
% these shocks are IID so this is true for both
% initial states
omega = (paramvec(3,:));
% omega = paramvec(3,:);

rho = (paramvec(5,:));
% sigma = 0.5;
sigma = (paramvec(4,:));

v = (paramvec(6, :));

% rho = 0.25;
% sigma = paramvec(5,:);
% rho = paramvec(4,:);

% size of technology shock
% this implies an asymptotic
% mean of 1 for the xi process
kappa = g / omega;
% kappa = 0.01;
% kappa = 0.01;

% number of theta states
% used to approximate continuous density
n_gridpoints = 80;

% reversed exponential growth per Dimitris
% top_grid = - normcdf(paramvec(7,:)) * 5;
theta0 = (paramvec(7,:));
growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;
theta_grid = (theta0).*((1 + growth_rate).^(1:(n_gridpoints)));
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

A_1(2:end, 2:end) = A_0(2:end, 2:end) * (1 - alpha);
A_1(2:end,2) = A_1(2:end,2) + alpha;

% transpose for use w/ Bianchi formulas
% VAR format

A_0 = A_0';
A_1 = A_1';


% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde = A_0(1:(end-1),1:(end-1));
A_1_tilde = A_1(1:(end-1),1:(end-1));

c_0_tilde = [0; A_0(2:end-1,end)];
c_1_tilde = [kappa; A_1(2:end-1,end)];

A_0_tilde(2:end,2:end) = A_0_tilde(2:end,2:end)- repmat(c_0_tilde(2:end,1),1,size(c_0_tilde,1)-1);
A_1_tilde(2:end,2:end) = A_1_tilde(2:end,2:end)- repmat(c_1_tilde(2:end,1),1,size(c_1_tilde,1)-1);


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
A_1_no_delta(2:end,2) = A_1_no_delta(2:end,2) + alpha;

% transpose for use w/ Bianchi formulas
% VAR format

A_0_no_delta = A_0_no_delta';
A_1_no_delta = A_1_no_delta';


% next, we will define subsetted matrices which omit the final column (this
% imposes the restriction that probabilities sum to 1)
A_0_tilde_no_delta = A_0_no_delta(1:(end-1),1:(end-1));
A_1_tilde_no_delta = A_1_no_delta(1:(end-1),1:(end-1));

c_0_tilde_no_delta = [0; A_0_no_delta(2:end-1,end)];
c_1_tilde_no_delta = [kappa; A_1_no_delta(2:end-1,end)];

A_0_tilde_no_delta(2:end,2:end) = A_0_tilde_no_delta(2:end,2:end)- repmat(c_0_tilde_no_delta(2:end,1),1,size(c_0_tilde_no_delta,1)-1);
A_1_tilde_no_delta(2:end,2:end) = A_1_tilde_no_delta(2:end,2:end)- repmat(c_1_tilde_no_delta(2:end,1),1,size(c_1_tilde_no_delta,1)-1);


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
mu = w * q;
mu2 = mu * mu';

wtilde = [w, zeros((n_coefs-1),2)];
qtilde = [q; piVec];

steady_state = [mu(2:(n_coefs-1)); 1 - sum(mu(2:end))];

figure(1)
plot(theta_grid,steady_state, "-o")
title('Steady State Theta Grid')


figure(9)
plot(log(theta_grid),steady_state, "-o")
title('Steady State Theta Grid, log scale for x axis')
xlabel('log(theta)')
% steady state values
H_star = theta_grid * steady_state;
L_star = 1 - H_star;
xi_star = mu(1);

% shock the system
% I'm 100% sure there's a cleaner way to do this but
% hey, it works
xi_shock = xi_star + kappa;
shock_state = steady_state - alpha * steady_state;
shock_state(1) = shock_state(1) + alpha;

shock_vec = [xi_shock;shock_state];
figure(2)
plot(shock_state)
xi_var = kappa^2 / (2 * g - g^2) * (1 - omega) * (omega);

% fsolve() for 2 x 2 system
% choose mu and lambda to match labor share = 0.66
% and ratio of high/low wages
options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
lambdamu = fsolve(@(x) calcloss(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta), [0; 0], options); 

% normcdf(lambdamu)

theor_mom = calcmom(lambdamu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v, ...
    A_0_tilde_no_delta, A_1_tilde_no_delta, c_0_tilde_no_delta, ...
    c_1_tilde_no_delta);

lambda = normcdf(lambdamu(1));
mu = normcdf(lambdamu(2));
% we target labor share above, now we target other xistuff
theor_mom = theor_mom(3:end,:);

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

emp_mom = [0.025; -0.0125; 0; 0; 0; ...
             emp_abs_wage_growth; ...
             emp_wage_growth; ...
             tenth_pctile_probs];
        
weight_vec = [5; 3; 1; 1; 1;
         1; 1; 1; 1; 1; ...
         3.5; 3; 3; 3; .5; ...
         2.5; 2; 2; 2; 2.5];

loss_vec = (theor_mom - emp_mom) ./ (0.01 + abs(emp_mom)) .* weight_vec;
loss_vec = [loss_vec; bottom_density_loss; top_density_loss];

% bars(labels, loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
% 
figure(2)
momlabels = categorical(1:17, 1:17, {'Output IRF','LShare IRF',...
     'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
     'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]',...
     'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]'},...
     'Ordinal',true);
bar(momlabels', [theor_mom([1:2, 6:end]), emp_mom([1:2, 6:end])])
title('Moment Matching (excluding signs)')

% 
labels = categorical(1:22, 1:22, {'Output IRF','LShare IRF',...
    'Output ','Wage Sign', 'Lshare IRF sign', ...
     'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
     'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]',...
     'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]',...
     'Bottom Density Penalty', 'Top Density Penalty'}, 'Ordinal',true);

figure(3)
bar(labels', loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
title('Weighted Percent Loss Contribution')


figure(4)
bar(labels(1:20)', [weight_vec .* weight_vec])
title('Moment Weights')

displaymat = [theor_mom, emp_mom];
% displaymat([1:2, 6:end], :)


   xi_shock = xi_star + kappa;
   shock_state = steady_state - alpha * steady_state;
   shock_state(1) = shock_state(1) + alpha;

   H_star = theta_grid * steady_state;
   L_star = 1 - H_star;
   
   
   % steady state production values
   X_star = calc_X(xi_star, L_star, lambda, rho);
   Y_star = calc_Y(H_star, X_star, mu, sigma, v);

   % high and low wages at steady state
   high_wage = w_h(H_star, L_star, xi_star, rho, sigma, mu, lambda, v);
   low_wage = w_l(H_star, L_star, xi_star, rho, sigma, mu, lambda, v);

   
   H = zeros(n_periods + 1, 1);
   L = zeros(n_periods + 1, 1);
   xi = zeros(n_periods + 1, 1);
   
   wh = zeros(n_periods + 1, 1);
   wl = zeros(n_periods + 1, 1);
   
   H(1) = H_star;
   L(1) = L_star;
   xi(1) = xi_star;
   wh(1) = high_wage;
   wl(1) = low_wage;
  
   H(2) = theta_grid * [shock_state];
   L(2) = 1 - H(2);
   xi(2) = xi_shock;
   wh(2) = w_h(H(2), L(2), xi(2), rho, sigma, mu, lambda, v);
   wl(2) = w_l(H(2), L(2), xi(2), rho, sigma, mu, lambda, v);
   
   
   shock_vec = [xi_shock; shock_state(1:(end - 1))];
   for i=1:(n_periods - 1)
       shock_vec = (1 - omega) * (A_0_tilde * shock_vec + c_0_tilde) + ...
           omega  * (A_1_tilde * shock_vec + c_1_tilde);
          
       xi(i + 2) = shock_vec(1,:);
       shock_state = shock_vec(2:end,:);

       shock_vec_for_calcs = [shock_state; 1 - sum(shock_state)];

       H(i + 2) = theta_grid * shock_vec_for_calcs;
       L(i + 2) = 1 - H(i + 2);

       % shock state production values
       X_shock = calc_X(xi(i + 2), L(i + 2), lambda, rho);
       Y_shock = calc_Y(H(i + 2), X_shock, mu, sigma, v);

       % high and low wages at shock state
       wh(i + 2) = w_h(H(i + 2), L(i + 2), xi(i + 2), rho, sigma, mu, lambda, v);
       wl(i + 2) = w_l(H(i + 2), L(i + 2), xi(i + 2), rho, sigma, mu, lambda, v);

   end
   
figure(5)
subplot(3,2,1);
plot([0, (1:20)./4]', (wh(1:end)./wh(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("High Wage")

subplot(3,2,2); 
plot([0, (1:20)./4]', (wl(1:end)./wl(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("Low Wage")

subplot(3,2,3);
plot([0, (1:20)./4]', (L(1:end)./L(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("L Skill Level")

subplot(3,2,4); 
plot([0, (1:20)./4]', (H(1:end)./H(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("H Skill Level")

subplot(3,2,5); 
plot([0, (1:20)./4]', (xi(1:end)./xi(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("Technology Level")

wage_diff = wh - wl;
subplot(3,2,6); 
plot([0, (1:20)./4]', (wage_diff(1:end)./wage_diff(1) - 1)/(kappa / sqrt(xi_var)), '.-')
title("High Wage - Low Wage")

names = {'HC Increase Prob', 'Conditional Fall Prob', ...
    'Shock Prob', 'Skilled Share', 'Unskilled Share', 'Skilled Curvature', ...
    'Unskilled Curvature', 'DRS Param', 'Bottom Rung'};

all_params = [phi, alpha, omega, mu, lambda, sigma, rho, v, theta0]';

disp(table(names', all_params))




