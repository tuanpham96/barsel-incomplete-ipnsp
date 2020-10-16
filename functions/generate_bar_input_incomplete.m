function X = generate_bar_input_incomplete(N, p, p_inc, normX)
X = zeros(N); 

index_on = rand(N,2) < p;

X(index_on(:,1),:) = 1;
X(:,index_on(:,2)) = 1;

X = X(:); 

index_on = find(X); 
n_incomplete = round(length(index_on) * p_inc); 
X(index_on(randperm(length(index_on), n_incomplete))) = 0;

if normX
    X_norm = sqrt(sum(X.^2));
    
    if X_norm > 0
        X = X / X_norm;
    end
end
end