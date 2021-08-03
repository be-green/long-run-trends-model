function end_grid = single_period_change(grid, s0_trans, s1_trans,...
    s0_intercept, s1_intercept, omega)
    short_grid = grid(1:(end - 1));
    moving_grid = (1 - omega) * (s0_trans * short_grid + s0_intercept) ...
        + omega * (s1_trans * short_grid + s1_intercept);
    end_grid = [moving_grid; 1 - sum(moving_grid)];
end