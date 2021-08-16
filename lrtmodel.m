function loss = lrtmodel(paramvec)

% transpose for particleswarm
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

% unconditional growth

% number of iterations associated with a single shock
% since a shock is 5 years, this corresponds to there being
% 4 iterations of the VAR a year
n_periods = 60;

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
sigma = rho - (paramvec(4,:));

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

if theta_order == 0
    % growth_rate = paramvec(8,:);
    % n_gridpoints = floor(-log(theta0) / log(1 + growth_rate));
    growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;
    theta_grid = (theta0).*((1 + growth_rate).^(1:(n_gridpoints)));
    
else 
    growth_rate = exp((log(theta0)) / n_gridpoints) - 1;
    theta_grid = 1 - (1.*((1 + growth_rate).^(1:(n_gridpoints))));
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
% figure(1)
% plot(theta_grid,steady_state)
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
% figure(2)
% plot(shock_state)
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
    c_1_tilde_no_delta, 1);

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
alphatomega = (alpha * omega < 0.0028) * abs(alpha * omega - 0.0028)/0.0028 * 100;
% half_income_from_low_skill = low_wage *(1-theta_0) / (low_wage * (1-theta_0) +  theta_0 * high_wage) >= 0.5;


emp_mom = [0.0281; -0.0125; 0; 0; 0; ...
             emp_abs_wage_growth; ...
             emp_wage_growth; ...
             tenth_pctile_probs];
        
weight_vec = [20; 10; 1; 1; 1;
         1; 1; 1; 1; 1; ...
         7; 4; 4; 4; 7; ...
         0; 0; 0; 0; 0];

loss_vec = (theor_mom - emp_mom) ./ (0.01 + abs(emp_mom)) .* weight_vec;
%  bars(labels, loss_vec .* loss_vec ./ (loss_vec' * loss_vec))
%  
%   figure(2)
%   momlabels = categorical(1:17, 1:17, {'Output IRF','LShare IRF',...
%          'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
%           'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]',...
%           'P(10)[0,25]','P(10)[25,50]','P(10)[50,75]','P(10)[75,95]','P(10)[95,100]'},...
%           'Ordinal',true);
%      bar(momlabels', [theor_mom([1:2, 6:(end - 5)]), emp_mom([1:2, 6:(end - 5)])])
%     title('Moment Matching (excluding signs)')
loss_vec = [loss_vec; bottom_density_loss; top_density_loss; alphatomega];

% miss = ([theor_mom([1:2, 6:end]) - emp_mom([1:2, 6:end])] ./ (0.01 + abs(emp_mom([1:2, 6:end])))).^2;
% miss = miss .* weight_vec ./ sum(miss .* weight_vec);
% 
% figure(3)
% bar(labels', miss)
% bar(labels', loss_vec .* loss_vec ./ (loss_vec' * loss_vec))


loss = loss_vec' * loss_vec;

