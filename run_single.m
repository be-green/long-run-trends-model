best_so_far = [0.183652720901365         0.405891342686913       0.00313580735609862         0.746562414814953         0.561574672485603       ...
    1.28866079486553         0.999926711104877         0.114394689025191         0.271698852278851          1.23898806841043 ...
    0.00737468986106264         0.374998460676595         0.183652720901365         0.920899983031557        0.0214624595047594];

other_good_one = [0.377582194256954          0.52280792265876       0.00424977207275125  ...
    0.323136129666296         0.372127259046202 0.740463682378398         0.999999666528352         0.329920455077772 ...
    0.538379273956443          1.49193350083659 0.0105572133006224         0.373385420269854         0.101850219127006   ...
    0.999999375290123        0.0833302279815142];

best_four_fixed_params = [0.184408030959959         0.358405991124413         0.003162510236958    ...
    0.910380774189953          0.63284847619654 ...
    0.083266759337691          0.29518975804057         0.271202225284933       ...
    0.00737468986106264         0.374998460676595 ...
    1        0.0214624595047594];
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
[csol, closs] = patternsearch(@(x) lrtmodel(x, 0, hyperparams), ...
                                        other_good_one, ...
                                        Aineq, bineq, [], [], ...
                                        lower, ...
                                        upper,...  
                                        [], patternoptions);
