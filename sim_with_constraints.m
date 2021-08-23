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
else
    error('Parsing function not yet adapted.')
end  

end