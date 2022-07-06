n_gridpoints = 120;
scale_period = 12;
n_periods = 1;
nstarts = 100;
parse_fcn_name = 'parse_model_params_v5';


weight_vec = [30; 0; 25; 25; 1; 1; 1;... labor share, wage ratio, labor share IRF, output IRF, % 3 sign restrictions
         0; 0; 0; 0; 0; ... abs wage moments
         15; 8; 8; 8; 30; ... wage moments
         40; ... wage difference between 5 and 4
         0; 0; 0; 0; 0; ... E(awg | income)
         0; 0; 0; 0; 0; ... E(wg | income)
         0; 0; ... E(awg), E(wg)
         5; 5; 5; 5; 5; ...
         6; ... % p10(5) - p10(1)
         0]; % aggregate standard deviation / sqrt(60)

hyperparams = struct('theta0', 0.03, 'scale_period', scale_period, ...
    'n_gridpoints', n_gridpoints, 'n_periods', n_periods, 'H_inside', 0, ...
    'parse_fcn_name', parse_fcn_name, 'weight_vec', weight_vec);

best_five_fixed_params = [0.184408030959959         0.358405991124413         0.003162510236958    ...
    0.910380774189953          0.63284847619654 ...
    0.083266759337691          0.29518975804057         0.271202225284933       ...
    0.00737468986106264        0.0214624595047594];
    
lrtmodel(best_five_fixed_params, 1, hyperparams)