timestamp = '23-Aug-2021 21_56_44';
function_vals = [];
sols = [];
exit_flags = [];

for i = 1:1500
   try
       % this should have worked but there was a bug
       % this wil be a lot faster after we fix it
%       load(['./model-output/model-run-number',num2str(i),'/runfeedback.mat']) 
%       function_vals = [function_vals; i, loss, exit];
%       sols = [sols; sol];
    
       close all     
       % so we can hack it by running publish code
        outdir = ['./model-output_',timestamp, '/model-run-number',num2str(i)];
%         addpath(outdir);
        
        load([outdir,'/runfeedback.mat'])
%         fid = fopen([outdir, '/publishcode.m']);
%         % we want the second line, this gets it
%         line = fgetl(fid);
%         line = fgetl(fid);
%         fclose(fid);
%         
%         eval(line);
%         loss = lrtmodel(this_solution, 0, 0, 80, 'parse_model_params_v2');
        
        function_vals = [function_vals; i, loss];
        sols = [sols; sol];
        exit_flags = [exit_flags;exit];
        rmpath(outdir);
        clear loss sol exit
   catch    
   end 
end    


% note that we can save this off in a folder if we want...
% basic idea would be to make some numbered folders called best1 -- best10
% then, we copy over the publish script there, then compile the pdf in each
% folder

trunc_fvals = function_vals;
trunc_sols = sols;

top_sols = [];
top_fvals = [];

% create a folder to store the best outputs
outdir = ['./model-output_', timestamp, '/best-solutions'];
addpath(outdir);

if ~exist(outdir, 'dir')
   mkdir(outdir)
end

for j = 1:5
    display(['Plotting/saving output for point rank number ',num2str(j)])

    best_f = min(trunc_fvals(:,2));
    best_sol = trunc_sols(trunc_fvals(:,2)==min(trunc_fvals(:,2)),:);
    
    best_i = trunc_fvals(trunc_fvals(:,2)==min(trunc_fvals(:,2)),1)
    
    best_index_vec = trunc_fvals(:,2)==min(trunc_fvals(:,2));
    
    top_fvals = [top_fvals; best_f];
    top_sols = [top_sols; best_sol];
    
    trunc_fvals(best_index_vec,:) = []; % delete best one;
    trunc_sols(best_index_vec,:) = []; % delete best one;
    
    publishpath = ['./model-output_', timestamp, '/model-run-number',num2str(best_i)];
    
        
    addpath(publishpath)
    if exist([publishpath, '/publishcode.pdf'])
        delete [publishpath, '/publishcode.pdf'];
    end
    theseoptions = struct('showCode', false, 'format','pdf','outputDir',publishpath,...
                           'codeToEvaluate','parse_model_params_v2 = ''parse_model_params_v2'';');

    publish([publishpath, '/publishcode.m'], theseoptions);
    rmpath(publishpath);
    
    % now, copy over the output file
     inputfile = [publishpath,'/publishcode.pdf'];

     outputfile = [outdir,'/publishcode_best',num2str(j),'.pdf'];
     copyfile(inputfile,outputfile);  

end    

close all
hist(function_vals(function_vals(:,2) < 1000,2),50)

% looking at dispersion in alpha across the solutions
scatter(function_vals(function_vals(:,2) < 500,2),sols(function_vals(:,2) < 500,2))

scatter(function_vals(function_vals(:,2) < 500,2),sols(function_vals(:,2) < 500,3))


close all

% regenerate output for best solution
n_gridpoints=80;
trial = top_sols(1,:);
trial(6) = trial(6)*10;
lrtmodel(trial, 0, 1, n_gridpoints, 'parse_model_params_v2');
