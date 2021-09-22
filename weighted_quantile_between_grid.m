function value = weighted_quantile_between_grid(x, w, prob)
    value = wprctile(x, prob * 100, w / sum(w));
end