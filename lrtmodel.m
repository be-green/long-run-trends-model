function loss = lrtmodel(paramvec)

% transpose for particlefilter
paramvec = paramvec'

% set up variables in accordance w/ our model

% param vec = [phi; alpha; omega; rho; sigma;]

% rate of moving up the ladder
phi = normcdf(paramvec(1,:));
% phi = paramvec(1,:);

% conditional fall probability
alpha = normcdf(paramvec(2,:));
% alpha = paramvec(2,:);

% unconditional growth
g = 0.02;

% number of iterations associated with a single shock
% since a shock is 5 years, this corresponds to there being
% 4 iterations of the VAR a year
n_periods = 20;

% death rate
% corresponds to an average 40 year working life
delta = 0.025;

% arrival rate of technology
% in this case this is also the probability
% at period t
% of state 1 happening in period t + 1
% these shocks are IID so this is true for both
% initial states
omega = normcdf(paramvec(3,:));
% omega = paramvec(3,:);

rho = 1 - exp(paramvec(5,:))
% sigma = 0.5;
sigma = rho - 0.1 - exp(paramvec(4,:))

v = normcdf(paramvec(6, :));

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
n_gridpoints = 200;

% reversed exponential growth per Dimitris
% top_grid = - normcdf(paramvec(7,:)) * 5;
top_grid = 0;
thetabar = 0;
theta_grid = flip(1 - exp(linspace(-5, top_grid, n_gridpoints)));
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
% plot(steady_state)
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
lambdamu = fsolve(@(x) calcloss(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v), [0; 0]); 

normcdf(lambdamu)

theor_mom = calcmom(lambdamu, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v);

% we target labor share above, now we target other xistuff
theor_mom = theor_mom(3:end,:);

emp_abs_wage_growth = ... 
    [-0.002551; ...
     0.0009579; ...
     0.004768; ...
     0.007638; ...
     0.02014];

emp_wage_growth = [-0.01486; -0.01008; -0.01178; -0.01167; -0.02467];

emp_mom = [0.025; -0.0125; 1; -1; ...
            -1;
             emp_abs_wage_growth; ...
                 emp_wage_growth];
        
weight_vec = [10; 4; 1; 1; 1;
         0; 0; 0; 0; 0; ...
         5; 2; 2; 2; 5];


loss_vec = (theor_mom - emp_mom) ./ (0.01 + abs(emp_mom)) .* weight_vec;
% 
figure(2)
labels = categorical(1:12, 1:12, {'Output IRF','LShare IRF',...
     'AWG[0,25]','AWG[25,50]','AWG[50,75]','AWG[75,95]','AWG[95,100]', ...
     'WG[0,25]','WG[25,50]','WG[50,75]','WG[75,95]','WG[95,100]'}, 'Ordinal',true);
bar(labels', [theor_mom([1:2, 6:end]), emp_mom([1:2, 6:end])])

miss = ([theor_mom([1:2, 6:end]) - emp_mom([1:2, 6:end])] ./ (0.01 + abs(emp_mom([1:2, 6:end])))).^2;
miss = miss ./ sum(miss);
% 
% figure(3)
% bar(labels', miss)
% [theor_mom([1:2, 6:end]), emp_mom([1:2, 6:end])]

displaymat = [theor_mom, emp_mom];
% displaymat([1:2, 6:end], :)

% paramvec

loss = loss_vec' * loss_vec;

