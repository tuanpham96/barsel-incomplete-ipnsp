function [distrib, pairs] = get_binary_pairwise_measure_distribution(X, measure_fun)
num_words = size(X, 1); 
pairs = nchoosek(1:num_words, 2);  
distrib = nan(size(pairs,1), 1);  
for ind = 1:length(distrib)
    pw = pairs(ind, :); 
    X1 = X(pw(1), :); 
    X2 = X(pw(2), :);
    distrib(ind) = measure_fun(X1, X2);
end
end