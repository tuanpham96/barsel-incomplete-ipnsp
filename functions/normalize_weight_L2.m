function W = normalize_weight_L2(W)
W = W ./ sqrt(sum(W.^2, 2)); 
end