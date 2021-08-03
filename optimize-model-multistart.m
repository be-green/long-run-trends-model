prob = createOptimProblem('fminunc', 'objective', @lrtmodel, ...
      'x0', [-1.6684 -0.5500 -1.0599 -1.0585 -0.7342 1.6934]');
%  
MS = MultiStart('UseParallel',1);
%  
[sol, loss] = run(MS, prob, 100);
%  
% lrtmodel([0.1, 0.05, 0.2, 0.1, 0.4, 0.4]')