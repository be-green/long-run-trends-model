function [ upper, lower, Aineq, bineq] = build_constraints_v1(n_periods,n_gridpoints)

g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;
% death rate
% corresponds to an average 40 year working life
delta =  exp(log(0.025 + 0.02 + 1) / (n_periods / 5)) - 1;
delta = delta + g; % is g being added twice here?, but maybe only used for the bounds
top_density_rho_const = 0.02^(1/n_gridpoints);

% 0 >= omega * alpha - alpha
bineq = 0;
Aineq = [0, -1, 1, 0, 0, 0, 0, 0, 0];

% top gridpoint has < 2% of the total mass
% implies top density rho constraint 
bineq = [bineq; -top_density_rho_const * delta];
Aineq = [Aineq; [(top_density_rho_const - 1), 0, top_density_rho_const, 0, 0, 0, 0, 0, 0]];

% bottom gridpoint has < 10% of the total mass
% implies botom density constraint , where
% -0.9 delta >= 0.9 alpha * omega - 0.1 * phi
bineq = [bineq; -0.9 * delta];
Aineq = [Aineq; [-0.1, 0, 0.9, 0, 0, 0, 0, 0, 0]];

lower = [0.001, ... phi
    0.005, ... alpha
    1/360, ... omega * alpha
    0.25, ... H diff from L (curvature exponent)
    0,   ... L curvature exponent
    0.01, ... percent human capital loss
    0.1, ... a (for p_z log odds)
    0.001, ... lambda
    0.001]; % mu
upper = [0.3, ... % phi
         0.7,  ... % alpha
         1/12, ...    % omega
         3, ...H diff from L (curvature exponent)
         0.75,    ...L curvature exponent
         0.99, ... percent human capital loss
         5, ... a (for p_z log odds)
         0.999, ... lamba
         0.999]; % mu