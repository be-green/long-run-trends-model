run_number = strrep(datestr(datetime), ':', '_');
addpath('/home/software/knitro/12.0.0');
addpath('/home/software/knitro/12.0.0/knitromatlab');
options = knitro_options('KN_PARAM_FEASTOL',1e-10,'KN_PARAM_FEASTOLABS',1e-10, ...
    'KN_PARAM_OPTTOL', 1e-10, 'KN_PARAM_XTOL_ITERS', 20, 'KN_PARAM_GRADOPT', 3,...
    'KN_PARAM_HESSOPT', 6);

n_gridpoints = 80;
n_periods = 60;
nstarts = 1000;

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

% TODO: could check that objective function doesn't error at starting vals
startvals = sim_with_constraints(nstarts, upper, lower, Aineq, bineq, parse_fcn_name);
    if ~exist(['./model-output_',run_number], 'dir')
       mkdir(['./model-output_',run_number])
    end
save(['model-output_',run_number,'/starting-values.mat'], 'startvals')


% evaluate objective at one, just as a sanity check;
lrtmodel(startvals(1,:), 0, 1, n_gridpoints,parse_fcn_name)

sol = zeros(nstarts, length(lower));
loss = zeros(nstarts, 1);
exitflg = zeros(nstarts, 1);

parfor i = 1:nstarts
   [this_solution, this_loss, this_exit] = knitromatlab(@(x) lrtmodel(x, 0, 0, n_gridpoints, parse_fcn_name), ...
                                    startvals(i,:), ...
                                    Aineq, bineq, [], [], ...
                                    lower, ...
                                    upper,...
                                    [], [], ...
                                    options);
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
    fprintf(fid, ['parse_fcn_name = [''', 'parse_model_params_v3','''];\n' ]);
    fprintf(fid, ['this_solution = [',num2str(this_solution),'];\n' ]);
    fprintf(fid, ['n_gridpoints = [',num2str(n_gridpoints),'];\n' ]);
    fprintf(fid, 'lrtmodel(this_solution, 0, 1, n_gridpoints, parse_fcn_name);');
    fclose(fid);
    
    addpath(outdir);
    % publish([outdir, '/publishcode.m'], theseoptions)

    parsave(['./model-output_',run_number, '/model-run-number' ,num2str(i),'/runfeedback.mat'],...
        this_solution, this_exit,this_loss);
    close all
end
