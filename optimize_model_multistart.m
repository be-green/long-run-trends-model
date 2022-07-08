
n_gridpoints = 120;
scale_period = 12;
n_periods = 1;
nstarts = 2;
parse_fcn_name = 'parse_model_params_v6';

weight_vec = [30; 20; 25; 25; 1; 1; 1;... labor share, wage ratio, labor share IRF, output IRF, % 3 sign restrictions
         0; 0; 0; 0; 0; ... abs wage moments
         15; 8; 8; 8; 30; ... wage moments
         25; ... wage difference between 5 and 4
         0; 0; 0; 0; 0; ... E(awg | income)
         0; 0; 0; 0; 0; ... E(wg | income)
         0; 0; ... E(awg), E(wg)
         5; 5; 5; 5; 5; ...
         6; ... % p10(5) - p10(1)
         0]; % aggregate standard deviation / sqrt(60)
     
run_number = strrep(datestr(datetime), ':', '_');

fminconoptions = optimoptions(@fmincon,'MaxIterations',1000, 'Display', ...
    'iter', 'FiniteDifferenceType', 'central', 'ScaleProblem', 'obj-and-constr', ...
    'HessianApproximation', 'bfgs', 'StepTolerance', 1e-10, ...
    'MaxFunctionEvaluations', 20000, 'OptimalityTolerance', 1e-15, ...
    'InitBarrierParam', 0.5);
addpath('C:/Program Files/Artelys/');
addpath('C:/Program Files/Artelys/Knitro 12.4.0/');
addpath('C:/Program Files/Artelys/Knitro 12.4.0/knitromatlab/');

koptions = 'knitro_options.opt';

method = "patternsearch";

patternoptions = optimoptions('patternsearch','Display','iter','PlotFcn',[], ...
    'MaxIterations',2000, 'MaxFunctionEvaluations', 20000);

% n_gridpoints is number of gridopints used on the 0-1 interval
% scale_period is used to represent scaling of "rfsim5" measure in
% regressions
% n_periods is the number of periods used w/ IRFs
% so 1 indicates "impact responses"
% nstarts is # of starts given to multistart
% hyperparams is misc hyperparameters:
% (1) theta0: level of H at bottom rung of ladder

hyperparams = struct('theta0', 0.03, 'scale_period', scale_period, ...
    'n_gridpoints', n_gridpoints, 'n_periods', n_periods, 'H_inside', 0, ...
    'parse_fcn_name', parse_fcn_name, 'weight_vec', weight_vec);

% scale_period * 5 is because within build_constraints the scale factor is
% divided by 5
if strcmp(parse_fcn_name,'parse_model_params_v1')
    [ upper, lower, Aineq, bineq] = build_constraints_v1(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v2')
    [ upper, lower, Aineq, bineq] = build_constraints_v2(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v3')
    [ upper, lower, Aineq, bineq] = build_constraints_v3(scale_period * 5,n_gridpoints, hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v4')
    [ upper, lower, Aineq, bineq] = build_constraints_v4(hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v5')
    [ upper, lower, Aineq, bineq] = build_constraints_v5(hyperparams);
elseif strcmp(parse_fcn_name,'parse_model_params_v6')
    [ upper, lower, Aineq, bineq] = build_constraints_v6(hyperparams);
else
    error('Parse function not coded yet!')
end     

% TODO: could check that objective function doesn't error at starting vals
startvals = sim_with_constraints(nstarts, upper, lower, Aineq, bineq, parse_fcn_name);
    if ~exist(['./model-output_',run_number], 'dir')
       mkdir(['./model-output_',run_number])
    end
save(['model-output_',run_number,'/starting-values.mat'], 'startvals')
save(['model-output_',run_number,'/hyperparams.mat'], 'hyperparams')

% evaluate objective at one, just as a sanity check;
lrtmodel(startvals(1,:), 1, hyperparams)

sol = zeros(nstarts, length(lower));
loss = zeros(nstarts, 1);
exitflg = zeros(nstarts, 1);

parfor i = 1:nstarts
    if method == "knitro" 
       [this_solution, this_loss, this_exit] = knitro_nlp(@(x) lrtmodel(x, 0, hyperparams), ...
                                        startvals(i,:), ...
                                        Aineq, bineq, [], [], ...
                                        lower, ...
                                        upper,...  
                                        [], [], struct('par_numthreads', 0), ...
                                        koptions); % set knitro options
    elseif method == "patternsearch"
        [this_solution, this_loss, this_exit] = patternsearch(@(x) lrtmodel(x, 0, hyperparams), ...
                                        startvals(i,:), ...
                                        Aineq, bineq, [], [], ...
                                        lower, ...
                                        upper,...  
                                        [], patternoptions); % set fmincon options
    else 
        [this_solution, this_loss, this_exit] = fmincon(@(x) lrtmodel(x, 0, hyperparams), ...
                                        startvals(i,:), ...
                                        Aineq, bineq, [], [], ...
                                        lower, ...
                                        upper,...  
                                        [], fminconoptions); % set fmincon options
    end
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
    fprintf(fid, 'lrtmodel(sol, 1, hyperparams);');
    fclose(fid);
    
    addpath(outdir);
    % publish([outdir, '/publishcode.m'], theseoptions)

    parsave(['./model-output_',run_number, '/model-run-number' ,num2str(i),'/runfeedback.mat'],...
        this_solution, this_exit,this_loss, hyperparams);
    close all
end
