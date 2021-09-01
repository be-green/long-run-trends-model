run_number = strrep(datestr(datetime), ':', '_');

custom_iter = optimoptions(@fmincon,'MaxIterations',1000, 'Display', ...
    'iter', 'FiniteDifferenceType', 'central', 'ScaleProblem', 'obj-and-constr', ...
    'HessianApproximation', 'bfgs', 'StepTolerance', 1e-13, ...
    'MaxFunctionEvaluations', 20000);

% n_gripoints is number of gridopints used on the 0-1 interval
% scale_period is used to represent scaling of "rfsim5" measure in
% regressions
% n_periods is the number of periods used w/ IRFs
% so 1 indicates "impact responses"
% nstarts is # of starts given to multistart
% hyperparams is misc hyperparameters:
% (1) theta0: level of H at bottom rung of ladder
n_gridpoints = 80;
scale_period = 12;
n_periods = 1;
nstarts = 1000;
hyperparams = struct('theta0', 0.01);

parse_fcn_name = 'parse_model_params_v4';
% scale_period * 5 is because within build_constraints the scale factor is
% divided by 5
if strcmp(parse_fcn_name,'parse_model_params_v1')
    [ upper, lower, Aineq, bineq] = build_constraints_v1(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v2')
    [ upper, lower, Aineq, bineq] = build_constraints_v2(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v3')
    [ upper, lower, Aineq, bineq] = build_constraints_v3(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v4')
    [ upper, lower, Aineq, bineq] = build_constraints_v4(scale_period * 5,n_gridpoints, hyperparams);
else
    error('Parse function not coded yet!')
end     

% TODO: could check that objective function doesn't error at starting vals
startvals = sim_with_constraints(nstarts, upper, lower, Aineq, bineq, parse_fcn_name);
    if ~exist(['./model-output_',run_number], 'dir')
       mkdir(['./model-output_',run_number])
    end
save(['model-output_',run_number,'/starting-values.mat'], 'startvals')


% evaluate objective at one, just as a sanity check;
lrtmodel(startvals(1,:), 0, 1, n_gridpoints,parse_fcn_name, n_periods, scale_period, hyperparams)

sol = zeros(nstarts, length(lower));
loss = zeros(nstarts, 1);
exitflg = zeros(nstarts, 1);

parfor i = 1:nstarts
   [this_solution, this_loss, this_exit] = fmincon(@(x) lrtmodel(x, 0, 0, n_gridpoints, ...
       parse_fcn_name, n_periods, scale_period, hyperparams), ...
                                    startvals(i,:), ...
                                    Aineq, bineq, [], [], ...
                                    lower, ...
                                    upper,...
                                    [], custom_iter);
    sol(i,:) = this_solution;
    loss(i,:) = this_loss;
    exitflg(i,:) = this_exit;
    
    outdir = [['./model-output_',run_number],'/model-run-number',num2str(i)];
    
    if ~exist(outdir, 'dir')
       mkdir(outdir)
    end
    
    theseoptions = struct('format','pdf','outputDir',outdir,...
        'showCode', false);
    
    
    fid = fopen([outdir, '/publishcode.m'], 'wt');
    fprintf(fid, ['parse_fcn_name = [''', parse_fcn_name,'''];\n' ]);
    fprintf(fid, ['load(','''./model-output_',run_number, '/model-run-number' ,num2str(i),'/runfeedback.mat'')\n']);
    fprintf(fid, ['n_gridpoints = [',num2str(n_gridpoints),'];\n' ]);
    fprintf(fid, ['n_periods = [',num2str(n_periods),'];\n' ]);
    fprintf(fid, ['scale_factor = [',num2str(scale_period),'];\n' ]);
    fprintf(fid, 'lrtmodel(sol, 0, 1, n_gridpoints, parse_fcn_name, n_periods, scale_period, hyperparams);');
    fclose(fid);
    
    addpath(outdir);
    % publish([outdir, '/publishcode.m'], theseoptions)

    parsave(['./model-output_',run_number, '/model-run-number' ,num2str(i),'/runfeedback.mat'],...
        this_solution, this_exit,this_loss, hyperparams);
    close all
end
