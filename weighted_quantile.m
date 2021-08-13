function value = weighted_quantile(x, w, prob)
    w = w./sum(w);
   [s, i] = sort(x);
   ws = w(i);
   
   total_mass = 0;
   index = 1;
   while total_mass < prob
       total_mass = total_mass + ws(index);
       index = index + 1;
   end
   
   value = s(index);
   
end