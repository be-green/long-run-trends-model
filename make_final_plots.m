n_gridpoints = 120;
scale_period = 12;
n_periods = 1;
nstarts = 100;
parse_fcn_name = 'parse_model_params_v5';

weight_vec = [30; 0; 25; 25; 1; 1; 1;... labor share, wage ratio, labor share IRF, output IRF, % 3 sign restrictions
         0; 0; 0; 0; 0; ... abs wage moments
         25; 8; 8; 20; 40; ... wage moments
         40; ... wage difference between 5 and 4
         0; 0; 0; 0; 0; ... E(awg | income)
         0; 0; 0; 0; 0; ... E(wg | income)
         0; 0; ... E(awg), E(wg)
         5; 5; 5; 5; 19; ...
         20; ... % p10(5) - p10(1)
         0]; % aggregate standard deviation / sqrt(60)

n_gridpoints = 120;
scale_period = 12;
n_periods = 1;
parse_fcn_name = 'parse_model_params_v5';

hyperparams = struct('theta0', 0.03, 'scale_period', scale_period, ...
    'n_gridpoints', n_gridpoints, 'n_periods', n_periods, 'H_inside', 0, ...
    'parse_fcn_name', parse_fcn_name, 'weight_vec', weight_vec);

final_cal = [0.195103100011429   0.315190716313622       0.00245505658480031      ...
    0.85765784328175          0.73646992423748         0.172774815978316 ...
    0.276635064952635         0.312074581673651       0.00932781486106264        0.0412370249448035];

lrtmodel(final_cal, 1, hyperparams)
