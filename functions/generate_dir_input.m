function X = generate_dir_input(N, p)
X = false(N); 

index_on = rand(N,2) < p;

X(index_on(:,1),:) = true;
X(:,index_on(:,2)) = true;

shift_on = rand(2*N-1,1) < p; 
shift_on = find(shift_on) - N;

idmat_shifted_1 = false(N); 
if ~isempty(shift_on)
    idmat_shifted_1 = arrayfun(@(s) diag(true(N-abs(s),1),s), shift_on, 'uni', 0);
    idmat_shifted_1 = mean(cat(3,idmat_shifted_1{:}),3) > 0;
end

shift_on = rand(2*N-1,1) < p; 
shift_on = find(shift_on) - N;

idmat_shifted_2 = false(N); 
if ~isempty(shift_on)
    idmat_shifted_2 = arrayfun(@(s) rot90(diag(true(N-abs(s),1),s)), shift_on, 'uni', 0);
    idmat_shifted_2 = mean(cat(3,idmat_shifted_2{:}),3) > 0;
end

X = X | idmat_shifted_1 | idmat_shifted_2; 

X = X(:); 

norm_X = sqrt(sum(X.^2, 'all'));

if norm_X > 0
    X = X / norm_X;
end

end