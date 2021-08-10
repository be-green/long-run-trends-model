% calibrate model
% param vec = [phi; alpha; omega; rho; sigma; DRS; theta0]

% prob = createOptimProblem('fminunc', 'objective', @lrtmodel, ...
%      'x0', [0.2088 0.369 -0.7  -0.1766 -0.0530]'.*100);
%  
% MS = MultiStart('UseParallel',1);
%  
% [sol, loss] = run(MS, prob, 1000);
%  
% lrtmodel([0.1, 0.05, 0.2, 0.1, 0.4, 0.4]')

% this gets a promising u-shape:
% [-1.1452; -1.2892; -0.9310; -1.68; -0.0961]
% init_sigma = log(0.5 - 0.4);
% init_rho = log(0.05);
% [sol, loss] = fminsearch(@lrtmodel, [-1.25; -1.4972; -1.4017; init_rho; init_sigma]);
% [sol, loss] = fminsearch(@lrtmodel, [norminv([0.0482 0.0672 0.0805]), -0.5, -2]');
% [sol, loss] = fminsearch(@lrtmodel, [0.4   -1.2   norminv(0.035)   -0.183   -0.79092    0.86]')
options = optimoptions('particleswarm',...
'Display', 'iter', 'UseParallel', true,'SwarmSize', 200);
[sol, loss] = particleswarm(@lrtmodel, 7, ...
    [0.001, 0.001, 0.001, -1, 0.25, 0.5, 0.001], ...
    [0.5, 0.95, 0.5, 0.25, 0.8, 0.9, 0.9], options);

% old particleswarm w/ 100 gridpoints
% sol = [0.4988    0.0095    0.0659    0.1557    0.8000    0.8545    0.0010]

paramvec = sol;
publish("plot_lrt_model.m", 'format','pdf', "showCode", false)