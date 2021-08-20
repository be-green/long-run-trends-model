custom_iter = optimoptions(@fmincon,'MaxIterations',500, 'Display', ...
    'iter', 'FiniteDifferenceType', 'central', 'ScaleProblem', 'obj-and-constr', ...
    'HessianApproximation', 'lbfgs');

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

lower = [0.001, ... phi * (p - ( 1 - p) )
    0.005, ... alpha
    1/360, ... omega * alpha
    0.25, ... H diff from L (curvature exponent)
    0,   ... L curvature exponent
    0.5, ... DRS parameter
    0.001,... theta0
    0.1, ... a (for p_z log odds)
    0.001, ... lambda
    0.001]; % mu
upper = [0.3, ... % phi * (p - (1 - p) )
         0.7,  ... % alpha
         1/120, ...    % omega
         3, ...H diff from L (curvature exponent)
         0.75,    ...L curvature exponent
         1, ... phi
         0.25, ... theta0
         10, ... a (for p_z log odds)
         0.999, ... lamba
         0.999]; % mu

nstarts = 1000;
startvals = sim_with_constraints(nstarts, upper, lower, A, b);

save('model-output/starting-values.mat', 'startvals')

sol = zeros(nstarts, length(lower));
loss = zeros(nstarts, 1);
exitflg = zeros(nstarts, 1);

parfor i = 1:nstarts
   [this_solution, this_loss, this_exit] = fmincon(@(x) lrtmodel(x, 0, 0, n_gridpoints), ...
                                    startvals(i,:), ...
                                    A, b, [], [], ...
                                    lower, ...
                                    upper,...
                                    [],custom_iter);
    sol(i,:) = this_solution;
    loss(i,:) = this_loss;
    exitflg(i,:) = this_exit;
    
    
    outdir = ['./model-output/model-run-number',num2str(i)];
    
    if ~exist(outdir, 'dir')
       mkdir(outdir)
    end
    
    theseoptions = struct('format','pdf','outputDir',outdir,...
        'showCode', false);
    
    
    fid = fopen([outdir, '/publishcode.m'], 'wt');
    fprintf(fid, ['this_solution = [',num2str(this_solution),'];\n' ]);
    fprintf(fid, ['n_gridpoints = [',num2str(n_gridpoints),'];\n' ]);
    fprintf(fid, 'lrtmodel(this_solution, 0, 1, n_gridpoints);');
    fclose(fid);
    
    addpath(outdir);
    publish([outdir, '/publishcode.m'], theseoptions)

    parsave(['./model-output/model-run-number',num2str(i),'/runfeedback.mat'],...
        this_solution, this_exit,this_loss);
    close all
end
