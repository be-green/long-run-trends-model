function end_wage_dist = single_bin_wage_dist_shock(i, C, shock_T, one_T, theta_grid)
    
    grid = zeros(size(theta_grid))';
    grid(i,:) = 1;
    
    % 100% chance of shock in first period
    grid = C + one_T * grid;
    grid = [grid; 1 - sum(grid)];

    
    % then iterate normally
    grid = C + shock_T * grid;
    grid = [grid; 1 - sum(grid)];
    
    end_wage_dist = grid;
end