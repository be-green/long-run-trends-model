function startvalues = sim_with_constraints(N, upper, lower, Aineq, bineq, parse_fcn_name)

if strcmp(parse_fcn_name,'parse_model_params_v1')

   which_unconstrained = sum(Aineq) == 0; 
   
   % weird bug means I can't do this directl
   sims = zeros(N, length(upper));
   for i = 1:N
       sims(i,:) = unifrnd(lower, upper).* which_unconstrained;
   end
   
   % sim relevant omega * alpha
   sims(:, 3) = unifrnd(lower(3), upper(3), [N, 1]);
   
   % divide by upper bound on omega
   alpha_lower = max(lower(2), sims(:,3) / 0.1);
   alpha_upper = upper(2);
   
   sims(:, 2) = unifrnd(alpha_lower, min(max(alpha_upper,alpha_lower+1e-4),1), [N, 1]);
   
   phi_upper = (0.9523 * sims(:,3) + bineq(1)) / 0.0477;
   phi_lower = (0.9 * sims(:,3) + bineq(2)) / 0.1;
   
   sims(:, 1) = unifrnd(phi_lower, phi_upper, [N, 1]);
   
   startvalues = sims;
elseif strcmp(parse_fcn_name,'parse_model_params_v2')
   which_unconstrained = min(abs(Aineq) == 0)+0; 
   
   % loop over obs so we can add a rejection sampling step
   sims = zeros(N, length(upper));
   for i = 1:N
       sims(i,:) = unifrnd(lower, upper).* which_unconstrained;

       solution_found = 0;
       while solution_found == 0
           % sim relevant omega * alpha * d -- expected HC loss per period
           sims(i, 3) = unifrnd(lower(3), upper(3));

           % first, we have bounds on phi to consider
           phi_l = max((abs(bineq(1)) + sims(i, 3))/ abs(Aineq(1,1)),lower(1));
           phi_u = min((abs(bineq(2)) + sims(i, 3))/ abs(Aineq(1,1)),upper(1));
           sims(i, 1) = unifrnd(phi_l, phi_u);

           % next, we need to pick levels of alpha which ensure that d
           % satisfies the bounds
           % constr 3 gives lower bound for d, so upper bound on alpha
           alpha_u = min(sims(i,3) / abs(Aineq(3,2)),upper(2)); 
           alpha_l = max(sims(i,3) / abs(Aineq(4,2)),lower(2));
           sims(i, 2) = unifrnd(alpha_l, alpha_u);
           
           % next, we draw a value of kappa
           sims(i,6) = unifrnd(lower(6), upper(6));
           % note: for some reason I did bounds in reverse order...
           xi_constant_l = max(abs(bineq(6)) - sims(i,6)*abs(Aineq(6,6)),lower(10));
           xi_constant_u = min(abs(bineq(5)) - sims(i,6)*abs(Aineq(5,6)),upper(10));
           sims(i, 10) = unifrnd(xi_constant_l, xi_constant_u);
           
           if max(isnan(sims(i,:))) == 0 && xi_constant_l <= xi_constant_u && ...
              phi_l <= phi_u && alpha_l <= alpha_u
                solution_found = 1;
           end
       end
       
   end
       disp('Fraction of simulated starting vals satisfying each constraint')
       mean(Aineq*sims'-bineq < 0,2)
   
   startvalues = sims;
   elseif strcmp(parse_fcn_name,'parse_model_params_v3')
   which_unconstrained = min(abs(Aineq) == 0)+0; 
   
   % loop over obs so we can add a rejection sampling step
   sims = zeros(N, length(upper));
   for i = 1:N
       sims(i,:) = unifrnd(lower, upper).* which_unconstrained;

       solution_found = 0;
       while solution_found == 0
           % sim relevant omega * alpha * d -- expected HC loss per period
           sims(i, 3) = unifrnd(lower(3), upper(3));

           % first, we have bounds on phi to consider
           phi_l = max((abs(bineq(1)) + sims(i, 3))/ abs(Aineq(1,1)),lower(1));
           phi_u = min((abs(bineq(2)) + sims(i, 3))/ abs(Aineq(1,1)),upper(1));
           sims(i, 1) = unifrnd(phi_l, phi_u);

           % next, we need to pick levels of alpha which ensure that d
           % satisfies the bounds
           % constr 3 gives lower bound for d, so upper bound on alpha
           alpha_u = min(sims(i,3) / abs(Aineq(3,2)),upper(2)); 
           alpha_l = max(sims(i,3) / abs(Aineq(4,2)),lower(2));
           sims(i, 2) = unifrnd(alpha_l, alpha_u);
           
           % next, we draw a value of kappa
           sims(i,6) = unifrnd(lower(6), upper(6));
           % note: for some reason I did bounds in reverse order...
           % xi_constant_l = max(abs(bineq(6)) - sims(i,6)*abs(Aineq(6,6)),lower(10));
           % xi_constant_u = min(abs(bineq(5)) - sims(i,6)*abs(Aineq(5,6)),upper(10));
           % sims(i, 10) = unifrnd(xi_constant_l, xi_constant_u);
           
           if max(isnan(sims(i,:))) == 0 && ...
              phi_l <= phi_u && alpha_l <= alpha_u
                solution_found = 1;
           end
       end
       
   end
       disp('Fraction of simulated starting vals satisfying each constraint')
       mean(Aineq*sims'-bineq < 0,2)
   
   startvalues = sims;
 elseif strcmp(parse_fcn_name,'parse_model_params_v4')
   which_unconstrained = min(abs(Aineq) == 0)+0; 
   
   % loop over obs so we can add a rejection sampling step
   sims = zeros(N, length(upper));
   for i = 1:N
       sims(i,:) = unifrnd(lower, upper).* which_unconstrained;

       solution_found = 0;
       while solution_found == 0
           % sim relevant omega * alpha * d -- expected HC loss per period
           sims(i, 3) = unifrnd(lower(3), upper(3));

           % first, we have bounds on phi to consider
           phi_l = max((abs(bineq(1)) + sims(i, 3))/ abs(Aineq(1,1)),lower(1));
           phi_u = min((abs(bineq(2)) + sims(i, 3))/ abs(Aineq(1,1)),upper(1));
           sims(i, 1) = unifrnd(phi_l, phi_u);

           % next, we need to pick levels of alpha which ensure that d
           % satisfies the bounds
           % constr 3 gives lower bound for d, so upper bound on alpha
           alpha_u = min(sims(i,3) / abs(Aineq(3,2)),upper(2)); 
           alpha_l = max(sims(i,3) / abs(Aineq(4,2)),lower(2));
           sims(i, 2) = unifrnd(alpha_l, alpha_u);
           
           % next, we draw a value of kappa
           sims(i,6) = unifrnd(lower(6), upper(6));
           % note: for some reason I did bounds in reverse order...
           % xi_constant_l = max(abs(bineq(6)) - sims(i,6)*abs(Aineq(6,6)),lower(10));
           % xi_constant_u = min(abs(bineq(5)) - sims(i,6)*abs(Aineq(5,6)),upper(10));
           % sims(i, 10) = unifrnd(xi_constant_l, xi_constant_u);
           
           if max(isnan(sims(i,:))) == 0 && ...
              phi_l <= phi_u && alpha_l <= alpha_u
                solution_found = 1;
           end
       end
       
   end
       disp('Fraction of simulated starting vals satisfying each constraint')
       mean(Aineq*sims'-bineq < 0,2)
   
   startvalues = sims;
else    
    error('Parsing function not yet adapted.')
end  

end