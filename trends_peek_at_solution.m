
function_vals = [];
sols = [];

for i = 1:1000
   try
       i
       % this should have worked but there was a bug
       % this wil be a lot faster after we fix it
%       load(['./model-output/model-run-number',num2str(i),'/runfeedback.mat']) 
%       function_vals = [function_vals; i, loss, exit];
%       sols = [sols; sol];
    
       close all     
       % so we can hack it by running publish code
        outdir = ['./model-output/model-run-number',num2str(i)];
        addpath(outdir);

        theseoptions = struct('showCode', false, 'format','pdf','outputDir',outdir);
        publish([outdir, '/publishcode.m'], theseoptions)
        
        function_vals = [function_vals; i, ans];
        sols = [sols; this_solution];
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
outdir = ['./model-output/best-solutions'];
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
    
    % now, copy over the output file
     inputfile = ['./model-output/model-run-number',...
                   num2str(best_i),'/publishcode.pdf'];

     outputfile = [outdir,'/publishcode_best',num2str(j),'.pdf'];
     copyfile(inputfile,outputfile);  

end    


close all

% regenerate output for best solution
lrtmodel(top_sols(1,:), 0, 1, n_gridpoints);
