function [pctile_id] = calc_single_bin_pctile_shock(i, C, shock_trans, ...
    one_trans, theta_grid, high_wage, low_wage, ...
    start_wage, pctile_wage)

    dist = single_bin_wage_dist_shock(i, C, shock_trans, one_trans, theta_grid);
    wage_grid = theta_grid * high_wage + (1 - theta_grid) * low_wage;
    wg = log(wage_grid') - log(start_wage);
    wg_lt_pctile = wg < pctile_wage;
    
    pctile_id = wg_lt_pctile' * dist;
end