function value = weighted_quantile(x, w, prob)
    w = w./sum(w);
   [s, i] = sort(x);
   ws = w(i);
   
   total_mass = 0;
   index = 1;
   total_mass = total_mass + ws(index);

   while total_mass < prob
       index = index + 1;
       total_mass = total_mass + ws(index);
   end
   
   value = s(index);
   
end