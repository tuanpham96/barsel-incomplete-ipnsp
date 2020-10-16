function [X, Y] = check_same_length(X, Y) 
X = X(:); Y = Y(:); 
if length(X) ~= length(Y), error('The 2 vectors need to have the same length!'); end
end