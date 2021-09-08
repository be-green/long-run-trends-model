
%        (1): phi, the rate at which skilled workers move up the grid
%        (2): alpha, the probability of being displaced, given an
%             innovation arrives
%        (3): expected_hc_loss, d * alpha x omega (arrival rate of new innovations)
%        (4): CES parameter on L - CES parameter on H
%        (5): CES parameter on L 
%        (6): Expectation of xi
%        (7): p_z, capturing higher probability of exposed
%             workers being displaced
%        (8): mu, parameter on labor in outer nest
%        (9): lambda, parameter on labor in inner nest
%       (10): kappa_share_of_xi_mean, defined as
%               kappa * omega / (kappa * omega + xi_constant)
%       (11): g, the depreciation rate of xi (NOT annualized)
%       (12): p0_share, share of people who never move up the ladder, stuck
%       at H = 0.
%       (13): gamma, which is (2p - 1) * phi, where p is the probability
%       that a move on the ladder is in the upward direction
%       (14): nu, (v), the DRS parameter
%       (15): alpha * omega, the shock probability parameter
lrtmodel([0.6183    0.4487    0.0019    0.5354    ...
    0.8443    0.2115   0.4477    0.5900 ...   
    0.5    1.5000    0.0054    0.0083    ...
    0.2146    1.0000    0.0348],...
    1, hyperparams);
