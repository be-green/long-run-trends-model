prob = createOptimProblem('fmincon', 'objective', @lrtmodel, ...
      'x0', [0.0990    0.0273    0.1000    0.2500    0.7500    0.9000    0.0225    1.3780]',...
  [],[],[],[],[0.005, 0.005, 0.005, 0.25, -1,   0.5, 0.005, -10], ...
    [0.1,   0.99,  0.1,    5, 0.75,    0.9, 0.25, 10]);
%  
MS = MultiStart('UseParallel',1);
%  
[sol, loss] = run(MS, prob, 500);
%  
% lrtmodel([0.1, 0.05, 0.2, 0.1, 0.4, 0.4]')