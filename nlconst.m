function [c, ceq] = nlconst(x)
    phi = x(1);
    omegatalpha = x(3);
    alpha = x(2);
    omega = omegatalpha / alpha;
    delta = 0.0053;
    rho = (phi / (delta + phi - delta *  phi + omega *  alpha));
    
    
    bottom_density = (delta + alpha * omega) /  ...
        ( delta + phi - delta *  phi + omega *  alpha );
    c = zeros(5, 1);
    c(1) = omegatalpha / alpha - 1;
    c(2) = bottom_density - 0.1;
    c(3) = 1 - (1 - rho^79) / (1 - rho) * bottom_density - 0.03;
    c(4) = (1 - rho^79) / (1 - rho) * bottom_density - 1;
    c(5) = 0.9 * omega * alpha - 0.1 * phi - 0.9 * delta;
    ceq = [];
end