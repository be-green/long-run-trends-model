function end_wage_dist = single_bin_wage_dist(i, s0_trans, s1_trans, ...
    s0_intercept, s1_intercept, omega, n_periods, theta_grid)
    
    grid = zeros(size(theta_grid))';
    grid(i,:) = 1;

    for i=1:n_periods
        grid = single_period_change(grid, s0_trans, s1_trans, ...
            s0_intercept, s1_intercept,  omega);
    end
    
    end_wage_dist = grid;
end