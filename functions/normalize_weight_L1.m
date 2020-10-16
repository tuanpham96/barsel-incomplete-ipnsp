function W = normalize_weight_L1(W)
W = W ./ sum(W, 2); 
end