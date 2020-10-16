function W = normalize_weight_gen(W,n)
% n: order
W = W ./ (sum(W.^n, 2) .^ (1/n)); 
end