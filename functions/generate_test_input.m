function [Xs, labels, p_inc_vec] = generate_test_input(num_data_per_pinc, N, p_sel, selrow, sigma_noise, p_incs, normX, shuffle_opt)
num_data = length(p_incs) * num_data_per_pinc; 
num_sel = round(num_data_per_pinc * p_sel); 

p_inc_vec = repmat(to_col_vec(p_incs), [1, num_data_per_pinc]);
labels = zeros(length(p_incs), num_data_per_pinc); 
labels(:,1:num_sel) = 1;

p_inc_vec = p_inc_vec(:)';
labels = labels(:)'; 

if shuffle_opt
    shuffled_ind = randperm(num_data);
    labels = labels(shuffled_ind);
    p_inc_vec = p_inc_vec(shuffled_ind);
end

Xs = arrayfun(@(label, p_inc) generate_1test(N, label, selrow, sigma_noise, p_inc, normX), ...
    labels, p_inc_vec, 'uni', 0);
Xs = horzcat(Xs{:}); 

end

function X = generate_1test(N, label, selrow, sigma_noise, p_inc, normX)

X = zeros(N);
n_complete = round(N * (1-p_inc));
X(selrow, randperm(N, n_complete)) = (label == 1);

X = X + sigma_noise * randn(size(X)); 
X = X(:); 

if normX
    X_norm = sqrt(sum(X.^2));
    
    if X_norm > 0
        X = X / X_norm;
    end
end


end