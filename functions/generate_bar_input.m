function X = generate_bar_input(N, p, normX)
X = zeros(N); 

index_on = rand(N,2) < p;

X(index_on(:,1),:) = 1;
X(:,index_on(:,2)) = 1;

X = X(:);

if normX
    X_norm = sqrt(sum(X.^2));
    
    if X_norm > 0
        X = X / X_norm;
    end
end

end