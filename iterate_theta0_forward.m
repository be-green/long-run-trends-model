n_gridpoints = 120;
scale_period = 12;
n_periods = 1;
parse_fcn_name = 'parse_model_params_v5';

best_so_far_fixed_param = [0.183652720901365         0.405891342686913       ...
    0.00313580735609862         0.746562414814953         0.561574672485603       ...
                0.114394689025191         ...
    0.271698852278851          0.222679276066183 ...
    0.00737468986106264         0.374998460676595         ...
             0.920899983031557        0.0214624595047594];

init_theta0 = 0.03;
x0 = best_so_far_fixed_param;
n_tries = 50;
theta0_sols = zeros(n_tries, size(best_so_far_fixed_param, 2));
losses = zeros(n_tries, 1);
theta0_list = zeros(n_tries, 1);

for i = 1:10
         theta0 = init_theta0 - 0.0005 * (i - 1);
hyperparams = struct('theta0', theta0, 'scale_period', scale_period, ...
    'n_gridpoints', n_gridpoints, 'n_periods', n_periods, 'H_inside', 0, ...
    'parse_fcn_name', parse_fcn_name, 'weight_vec', weight_vec);
[ upper, lower, Aineq, bineq] = build_constraints_v5(hyperparams);



[x0, closs] = patternsearch(@(x) lrtmodel(x, 0, hyperparams), ...
                                        reduced_test, ...
                                        Aineq, bineq, [], [], ...
                                        lower, ...
                                        upper,...  
                                        [], patternoptions);
losses(i,:) = closs;
theta0_sols(i, :) = x0;
theta0_list(i,:) = theta0;
end

best_sol = theta0_sols(losses == min(losses), :);
best_theta0 = theta0_list(losses == min(losses),:);
    
