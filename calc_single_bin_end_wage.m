function [wage_growth, abs_wage_growth] = calc_single_bin_end_wage(i, s0_trans, ...
    s1_trans,s0_intercept, s1_intercept, omega, theta_grid, n_periods, high_wage, low_wage, ...
    start_wage)

    dist = single_bin_wage_dist(i, s0_trans, s1_trans,...
        s0_intercept, s1_intercept, omega, n_periods, theta_grid);
    wage_grid = theta_grid * high_wage + (1 - theta_grid) * low_wage;
    wg = log(wage_grid') - log(start_wage);
    awg = abs(log(wage_grid') - log(start_wage));
    
    wage_growth = wg' * dist;
    abs_wage_growth = awg' * dist;
end