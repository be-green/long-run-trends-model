function [sol, loss] = coordinate_descent(init_vec, Aineq, bineq, lower, upper) 
    
    attempt = init_vec;
    maxiter = 100;
    num_params = size(init_vec, 2);
    iter = 0;
    while iter < maxiter
        iter = iter + 1;
        slot = mod(iter, 15);
        
         if (all(Aineq(:, slot) == 0))
            l = lower(:, slot); 
            u = upper(:, slot);
         else
             init_l = lower(:, slot);
             init_u = upper(:, slot);
             
             keep = repelem(1, num_params);
             keep(:, slot) = 0;
             Aineq_m_slot = Aineq(:, find(keep));
             attempt_m_slot = attempt(:, find(keep));
             Aineq_m_slot * attempt_m_slot' ./ 
                 
             
         end
        [attempt, thisloss] = patternsearch(@(x) run_search(x, slot, guess, hyperparams), [], [], [], [], [], l, u, patternoptions);
    
     end

    sol = attempt;
    loss = thisloss
end

function loss = run_search(guess, init_vec, slot, hyperparams)

    init_vec(:, slot) = guess;
    loss = lrtmodel(init_vec, 0, hyperparams);

end
        
        
    