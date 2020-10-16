function mi = calculate_mutual_information_discrete(X, Y, bin_size)
% Calculate mutual information of discrete distributions between X, Y
% (state vectors) as I(X,Y) = H(X) + H(Y) - H(X,Y) 

[X, Y] = check_same_length(X, Y); 
if length(bin_size) == 1, bin_size = [bin_size, bin_size]; end
if length(bin_size) ~= 2, error('"bin_size" has to either be a scalar or 2-element vector'); end

[~, edge_X] = bin_centers_and_edges(X,bin_size(1));
[~, edge_Y] = bin_centers_and_edges(X,bin_size(2));

entX = calculate_entropy_discrete(histcounts(X, edge_X, 'Normalization', 'probability'));
entY = calculate_entropy_discrete(histcounts(Y, edge_Y, 'Normalization', 'probability'));
entXY = calculate_entropy_discrete(histcounts2(X, Y, edge_X, edge_Y, 'Normalization', 'probability'));

mi = entX + entY - entXY;
end
