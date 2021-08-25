addpath('/home/software/knitro/12.0.0')
addpath('/home/software/knitro/12.0.0/knitromatlab')
custom_iter = optimoptions(@fmincon,'MaxIterations',500, 'Display', ...
    'iter', 'FiniteDifferenceType', 'central', 'ScaleProblem', 'obj-and-constr', ...
    'HessianApproximation', 'lbfgs');

n_gridpoints = 80;
n_periods = 60;

parse_fcn_name = 'parse_model_params_v3';

if strcmp(parse_fcn_name,'parse_model_params_v1')
    [ upper, lower, Aineq, bineq] = build_constraints_v1(n_periods,n_gridpoints);
elseif strcmp(parse_fcn_name,'parse_model_params_v2')
    [ upper, lower, Aineq, bineq] = build_constraints_v2(n_periods,n_gridpoints);
elseif strcmp(parse_fcn_name,'parse_model_params_v3')
    [ upper, lower, Aineq, bineq] = build_constraints_v3(n_periods,n_gridpoints);
else
    error('Parse function not coded yet!')
end     

[this_solution, this_loss, this_exit] = knitromatlab(@(x) lrtmodel(x, 0, 0, n_gridpoints, parse_fcn_name), ...
                                    [0.0820    0.1017    0.0029    0.8629    0.7529   10.6947    0.7564    0.2018    0.0521    0.0000    0.2000], ...
                                    Aineq, bineq, [], [], ...
                                    lower, ...
                                    upper,...
                                    []);