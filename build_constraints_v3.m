function [ upper, lower, Aineq, bineq] = build_constraints_v3(n_periods,n_gridpoints)

% g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;
% death rate
% corresponds to an average 40 year working life
delta =  exp(log(0.025 + 1) / (n_periods / 5)) - 1;
%delta = delta + g; % is g being added twice here?, but maybe only used for the bounds

omega = 5/n_periods; %fixing this ex ante in this round

% setting up theta grid
theta0 = 0.02;
growth_rate = exp((-log(theta0)) / n_gridpoints) - 1;

% set lower bound on expected human capital growth
dlogtheta_lb = 0.1 * delta; % expected level of HC at end of career is 0.1;

% set lower bound on expected human capital growth
dlogtheta_ub = 2 * delta; % expected level of HC at end of career is 0.6;

% set upper/lower bound on the implied human capital loss
min_d = 0.05;
max_d = 0.6;

% phi -  d* omega * alpha >= dlogtheta_lb
% ie -dlogtheta_lb <= d*omega*alpha - phi * growth_rate
bineq = -dlogtheta_lb;
Aineq = [-growth_rate, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% phi*growth_rate -  d* omega * alpha <= dlogtheta_ub
bineq = [bineq;dlogtheta_ub];
Aineq = [Aineq; growth_rate, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% d = paramvec(3) / (alpha * omega). So if we want d >= min_d, that is isomorphic to
% having alpha * omega * min_d <= paramvec(3), so -paramvec(3) +alpha * omega * min_d <= 0
bineq = [bineq;0];
Aineq = [Aineq; 0, min_d * omega, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0];

% d = paramvec(3) / (alpha * omega). So if we want d <= max_d, that is isomorphic to
% having (alpha * omega) * max_d >= paramvec(3), so -(alpha * omega) * min_d + paramvec(3) <= 0
bineq = [bineq;0];
Aineq = [Aineq; 0, -max_d * omega, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];


% set bounds on kappa to restrict that the mean of xi isn't insanely high...
% xi_mean = (omega kappa + xi_constant) / g
xi_mean_lb = 0.2;
xi_mean_ub = 5;

% only do upper bound for now?
% bineq = [bineq;xi_mean_ub*g];
% Aineq = [Aineq; 0, 0, 0, 0, 0, omega, 0, 0, 0, 1, 0];

% Adding the lower bound. Only do upper bound for now?
% bineq = [bineq;-xi_mean_lb*g];
% Aineq = [Aineq; 0, 0, 0, 0, 0, -omega, 0, 0, 0, -1, 0];

% standard deviation is also linear in kappa. Can impose that it's not huge
% relative to the mean of xi. But won't do this for now...
% xi_std = sqrt(kappa^2 / (2 * g - g^2) * (1 - omega) * (omega));

g_upper = 0.2;
g_lower = 0.02;

lower = [0.001, ... phi
    0.1, ... alpha
    1/12*0.005*min_d, ... d * omega * alpha (omega = 1 per year)
    0.25, ... H diff from L (curvature exponent)
    0,   ... L curvature exponent
    xi_mean_lb / omega * g_upper, ... size of xi jump kappa
    0.1, ... a (for p_z log odds)
    0.001, ... lambda
    0.001,... % mu
    0, ...% xi intercept
    g_lower, 0]; % g

upper = [0.4, ... % phi
         0.7,  ... % alpha
         1/12*0.7*max_d/4, ...    % d * omega * alpha (omega = 1 per year)
         3, ...H diff from L (curvature exponent)
         0.85,    ...L curvature exponent
         (xi_mean_ub -0.1) / omega * g_upper, ... size of xi jump
         5, ... a (for p_z log odds)
         0.999, ... lamba
         0.999,....% mu
         xi_mean_ub*g_upper / 2, ... % xi intercept
         g_upper, 0.1]; % g
     