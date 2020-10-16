function s = calculate_SMD(X,Y)
[X, Y] = check_same_length(X, Y); 
s = sum(xor(X, Y))/length(X); 
end
