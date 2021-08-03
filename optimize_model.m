% calibrate model
% param vec = [phi; alpha; omega; rho; sigma; kappa]

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
init_sigma = log(0.5 - 0.4);
init_rho = log(0.05);
% [sol, loss] = fminsearch(@lrtmodel, [-1.25; -1.4972; -1.4017; init_rho; init_sigma]);
% [sol, loss] = fminsearch(@lrtmodel, [norminv([0.0482 0.0672 0.0805]), -0.5, -2]');
% [sol, loss] = fminsearch(@lrtmodel, [0.4   -1.2   norminv(0.035)   -0.183   -0.79092    0.86]')
[sol, loss] = particleswarm(@lrtmodel, 6, ...
    [-3, -3, -3, -3, -3, -3], ...
    [3, 3, 3, 0.5, 1, 3]) 
