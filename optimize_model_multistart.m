custom_iter = optimoptions(@fmincon,'MaxIterations',200, 'Display', ...
    'iter', 'FiniteDifferenceType', 'central', 'ScaleProblem', 'obj-and-constr', ...
    'HessianApproximation', 'lbfgs');
% parpool()

n_gridpoints = 80;
n_periods = 60;

g = exp(log(0.02 + 1) / (n_periods / 5)) - 1;
% death rate
% corresponds to an average 40 year working life
delta =  exp(log(0.025 + 0.02 + 1) / (n_periods / 5)) - 1;
delta = delta + g;
top_density_rho_const = 0.02^(1/n_gridpoints);

% 0 >= omega * alpha - alpha
b = 0;
A = [0, -1, 1, 0, 0, 0, 0, 0, 0, 0];

% top gridpoint has < 2% of the total mass
% implies top density rho constraint 
b = [b; -top_density_rho_const * delta];
A = [A; [(top_density_rho_const - 1), 0, top_density_rho_const, 0, 0, 0, 0, 0, 0, 0]];

% bottom gridpoint has < 10% of the total mass
% implies botom density constraint , where
% -0.9 delta >= 0.9 alpha * omega - 0.1 * phi
b = [b; -0.9 * delta];
A = [A; [-0.1, 0, 0.9, 0, 0, 0, 0, 0, 0, 0]];

lower = [0.001, 0.005, 1/360, 0.25, -2,   0.5, 0.001, -10, 0.1, 0.1];
upper = [0.05,   0.7,  1/120,    5, 0.75,    1, 0.1, 10, 0.9, 0.9];

nstarts = 10;
startvals = sim_with_constraints(nstarts, upper, lower, A, b);

sol = zeros(nstarts, length(lower));
loss = zeros(nstarts, 1);
exitflg = zeros(nstarts, 1);

parfor i = 1:nstarts
   [sol(i,:), loss(i,:), exitflg(i,:)] = fmincon(@(x) lrtmodel(x, 0, 0, n_gridpoints), ...
                                    startvals(i,:), ...
                                    A, b, [], [], ...
                                    lower, ...
                                    upper,...
                                    [],custom_iter);
                                
                                publish(@(x) lrtmodel(x, 0, 1, n_gridpoints))
end


