function end_grid = single_period_change(grid, s0_trans, s1_trans,...
    s0_intercept, s1_intercept, omega)
    newgrid = (1 - omega) * (s0_trans * grid + s0_intercept) ...
        + omega * (s1_trans * grid + s1_intercept);
    % needs to sum to 1
    assert(sum(newgrid) - 1 < 1e-10, 'Probability density on grid failed to sum to 1 following transition.' );
    end_grid = newgrid;
end