custom_iter = optimoptions(@fmincon,'MaxIterations',200);

prob = createOptimProblem('fmincon', 'objective', @(x) lrtmodel(x, 0, 0, 80), ...
      'x0', [0.0990    0.0273    0.1000    0.2500    0.7500    0.9000    0.0225    1.3780],...
    'lb', [0.005, 0.005, 0.005, 0.25, -1,   0.5, 0.005, -10], ...
    'ub', [0.1,   0.99,  0.1,    5, 0.75,    1, 0.25, 10],...
    'options', custom_iter);
%  
MS = MultiStart('UseParallel',1, 'Display', 'iter', 'StartPointsToRun', 'bounds');
%  
[sol, loss] = run(MS, prob, 500);
%  
% lrtmodel([0.1, 0.05, 0.2, 0.1, 0.4, 0.4]')

function const = nlconst(alpha, omega, phi)
    const = (1 + phi) / (omega*alpha + 0.0053);
end
