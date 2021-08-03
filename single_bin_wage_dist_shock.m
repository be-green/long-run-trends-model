function end_wage_dist = single_bin_wage_dist_shock(i, s0_trans, s1_trans, ...
    s0_intercept, s1_intercept, omega, n_periods, theta_grid)
    
    grid = zeros(size(theta_grid))';
    grid(i,:) = 1;
    
    % 100% chance of shock in first period
    grid = single_period_change(grid, s0_trans, s1_trans, ...
        s0_intercept, s1_intercept, 1);
    
    % then iterate normally
    for i=1:(n_periods-1)
        grid = single_period_change(grid, s0_trans, s1_trans, ...
            s0_intercept, s1_intercept, omega);
    end
    
    end_wage_dist = grid;
end