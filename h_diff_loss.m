function loss = h_diff_loss(fudge_factor, omega_shock_amount)

final_cal = [0.195103100011429   0.315190716313622       0.00245505658480031      ...
    0.85765784328175          0.73646992423748         0.172774815978316 ...
    0.276635064952635         0.312074581673651       0.00932781486106264        0.0412370249448035];

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

% calculating this stuff to scale shocks correctly


% conditional fall probability
alpha = (final_cal(:, 2));

% arrival rate of technology
% in this case this is also the probability
% at period t
% of state 1 happening in period t + 1
% these shocks are IID so this is true for both
% initial states

% use alpha to get omega from alpha_x_omega
alpha_x_omega = final_cal(:, 10);
omega = alpha_x_omega / alpha;

% get d from d_omega_alpha / omega_alpha
d_x_omega_x_alpha = (final_cal(:, 3));
d = d_x_omega_x_alpha / (alpha * omega);
growth_rate = exp((-log(hyperparams.theta0)) / n_gridpoints) - 1;

% this is an increase of phi to balance the omega shock (keep H constant)
% we can approximate this w/ the process theta (g * phi - alpha omega d)
% if omega goes up by (omega + x), we need phi to increase by xd / g
% this approximation is a touch too low though, adding a fudge factor
% so as to match the steady state H from the the original calibration
phi_shock_amount = omega_shock_amount * alpha * d / growth_rate + fudge_factor;
[both_shock_loss, both_shock_emp_mom, both_shock_theor_mom, ...
    both_shock_wh, both_shock_wl, both_shock_Y, both_shock_H, both_shock_lshare,...
    both_shock_xi, both_shock_twenty_fifth_pctile, both_shock_seventy_fifth_pctile] = ...
    lrtmodel_shock_omega_phi(final_cal, 0, hyperparams, omega_shock_amount, phi_shock_amount);

loss = abs(both_shock_H(1) - both_shock_H(end)) * 1e10;
end


