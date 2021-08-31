function [ upper, lower, Aineq, bineq] = build_constraints_v4(n_periods,n_gridpoints, hyperparams)

% g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;
% death rate
% corresponds to an average 40 year working life
delta =  exp(log(0.025 + 1) / (n_periods / 5)) - 1;
%delta = delta + g; % is g being added twice here?, but maybe only used for the bounds

omega = 5/n_periods; %fixing this ex ante in this round

% setting up theta grid
theta0 = hyperparams.theta0;
growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;

% set lower bound on expected human capital growth
dlogtheta_lb = 0.1 * delta; % expected level of HC at end of career is 0.1;

% set lower bound on expected human capital growth
dlogtheta_ub = 2 * delta; % expected level of HC at end of career is 0.6;

% set upper/lower bound on the implied human capital loss
min_d = 0.05;
max_d = 0.6;

% growth_rate * phi * (2p_up - 1) -  d* omega * alpha >= dlogtheta_lb
% ie -dlogtheta_lb <= d*omega*alpha - phi * growth_rate
bineq = -dlogtheta_lb;
Aineq = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -log(1 + growth_rate)];

% growth_rate * phi * (2p_up - 1) -  d* omega * alpha <= dlogtheta_ub
bineq = [bineq;dlogtheta_ub];
Aineq = [Aineq; 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, log(1 + growth_rate)];

% d = paramvec(3) / (alpha * omega). So if we want d >= min_d, that is isomorphic to
% having alpha * omega * min_d <= paramvec(3), so -paramvec(3) +alpha * omega * min_d <= 0
bineq = [bineq;0];
Aineq = [Aineq; 0, min_d * omega, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% d = paramvec(3) / (alpha * omega). So if we want d <= max_d, that is isomorphic to
% having (alpha * omega) * max_d >= paramvec(3), so -(alpha * omega) * min_d + paramvec(3) <= 0
bineq = [bineq;0];
Aineq = [Aineq; 0, -max_d * omega, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% gamma = phi * (2 p_up - 1);
% we want 2 * p_up - 1 < 1, which means 
% gamma - phi < 0
% the other constraint (2 * p_up - 1) > 0 is given by gamma > 0, so 
% this is just a bound constraint
bineq = [bineq;0];
Aineq = [Aineq; 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];


lower = [0.001, ... phi
    0.1, ... alpha
    1/12*0.005*min_d, ... d * omega * alpha (omega = 1 per year)
    0.25, ... H diff from L (curvature exponent)
    0,   ... L curvature exponent
    0.2 , ... Mean of xi
    0.1, ... a (for p_z log odds)
    0.001, ... lambda
    0.001,... % mu
    0.5, ...% fraction of xi mean coming from kappa shocks vs intercept
    0.02/12, ... % g
    0 ... % p0 share
    0]; % p_up, conditional prob direction on ladder is up given move

upper = [0.4, ... % phi
         0.7,  ... % alpha
         1/12*0.7*max_d/4, ...    % d * omega * alpha (omega = 1 per year)
         3, ...H diff from L (curvature exponent)
         0.85,    ...L curvature exponent
         5, ... Mean of xi
         5, ... a (for p_z log odds)
         0.999, ... lamba
         0.999,....% mu
         1, ... % fraction of xi mean coming from kappa shocks vs intercept
         0.2/12, ... g
         0.25 ...  % p0 share
         1]; % p_up, conditional prob direction on ladder is up given move
     