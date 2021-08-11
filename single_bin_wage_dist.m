function end_wage_dist = single_bin_wage_dist(i, C, ss_trans, theta_grid)
    
    grid = zeros(size(theta_grid))';
    grid(i,:) = 1;
    
    grid = C + ss_trans * grid;
    grid = [grid; 1 - sum(grid)];
    
    end_wage_dist = grid;
end