
function loss = calcloss(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v)

   mom = calcmom(x, theta_grid, steady_state, xi_star, ...
    kappa, rho, sigma, alpha, phi, xi_var, A_0_tilde, A_1_tilde, ...
    c_0_tilde, c_1_tilde, omega, n_periods, v);
   % 12.2 taken via top 5% / bottom 5% regular wages
   theor_vec = mom(1:2,:);
   emp_vec = [0.66; 2.45];
   loss_vec = (theor_vec - emp_vec) ./ abs(emp_vec);
   loss = loss_vec' * loss_vec;
end