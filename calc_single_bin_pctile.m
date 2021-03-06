function [pctile_id] = calc_single_bin_pctile(i, s0_trans, ...
    s1_trans,s0_intercept, s1_intercept, omega, theta_grid, n_periods, high_wage, low_wage, ...
    start_wage, pctile_wage)

    dist = single_bin_wage_dist(i, s0_trans, s1_trans,...
        s0_intercept, s1_intercept, omega, n_periods, theta_grid);
    wage_grid = theta_grid * high_wage + (1 - theta_grid) * low_wage;
    wg = log(wage_grid') - log(start_wage);
    wg_lt_pctile = wg < pctile_wage;
    
    pctile_id = wg_lt_pctile' * dist;
end