function res = train_with_incomplete(opts) 

%% Network props 
N = opts.N; 
N_U = N^2;

a = opts.a_init; 
b = opts.b_init; 

mu = opts.mu;

eta_ip_a = opts.eta_ip_a;
eta_ip_b = opts.eta_ip_b;
eta_sp_ltp = opts.eta_sp_ltp; 
eta_sp_ltd = opts.eta_sp_ltp * opts.kappa_sp;

norm_W_order = opts.weight_properties.norm_order;
allow_only_exc = opts.weight_properties.onlyexcitatory; 

W = zeros(N); 
W(opts.selrow,:) = 1; 
W = W(:)' + rand(1,N_U); 

if norm_W_order > 0
    W = normalize_weight_gen(W, norm_W_order);
end 

reshaped_W = reshape(W,[N,N]);
bar_alphaW_0 = sum(reshaped_W(opts.selrow,:),'all') / sum(W(:));

%% Training generation
num_train = opts.num_train; 
p_train_complete = opts.p_train_complete;
num_train_complete = ceil(p_train_complete * num_train);
num_train_incomplete = num_train - num_train_complete; 
train_norm_input_opt = opts.train_norm_input;

p_bar = 1/(2*N); 
p_inc = opts.p_inc; 

vec_train_completeness = [ones(num_train_complete,1); (1-p_inc) * ones(num_train_incomplete,1)]'; 

U_trains = [arrayfun(@(~) generate_bar_input(N, p_bar, train_norm_input_opt),1:num_train_complete, 'uni', 0), ...
            arrayfun(@(~) generate_bar_input_incomplete(N, p_bar, p_inc, train_norm_input_opt), 1:num_train_incomplete, 'uni', 0)];

points_to_record_W = [1, num_train_complete, num_train]; 

%% Generate test data
test_p_inc_opts =  opts.test_p_incs;
[U_tests, test_labels, test_p_inc_vec] = generate_test_input(...
    opts.num_test_per_pinc, N, ...
    opts.test_p_sel, opts.selrow, ...
    opts.test_sigma_noise, test_p_inc_opts,...
    opts.test_norm_input, opts.test_shuffle);

%% Train and record 
subsampled = 1; 
if isfield(opts, 'subsampled')
    subsampled = opts.subsampled;
end

num_savedsamples = round(num_train/subsampled);

a_v = zeros(1, num_savedsamples);
b_v = zeros(1, num_savedsamples);
t_v = zeros(1, num_savedsamples);

W_s = cell(length(points_to_record_W), 1);

bar_alphaW_v = zeros(1, num_savedsamples);
dY_objvsnoise_test = zeros(length(test_p_inc_opts), num_savedsamples);

cnt_save = 1;
for i = 1:num_train
    % calculate from train
    U = U_trains{i}; 
    X = W * U; 
    Y = sigmoidal_activation(X, a, b);        
    
    % plasticity 
    [a, b] =  IP_Triesch_diffrate(X, Y, a, b, eta_ip_a, eta_ip_b, mu);
    W = SP_HebbMod(U, Y, W, eta_sp_ltp, eta_sp_ltd);
    
    if allow_only_exc
        W = max(W, 0);
    end
    
    if norm_W_order > 0
        W = normalize_weight_gen(W, norm_W_order);
    end
    
    % saving         
    reshaped_W = reshape(W,[N,N]); 
    if ~isempty(find(i == points_to_record_W,1))
        W_s{i} = reshaped_W;
    end
    
    if mod(i-1,subsampled) ~= 0 
        continue;
    end
    
    % testing
    X_tests = W * U_tests;
    Y_tests = sigmoidal_activation(X_tests, a, b);
    
    dY_objvsnoise_test(:,cnt_save) = arrayfun(@(p_inc_test) ...
        mean(Y_tests(test_p_inc_vec == p_inc_test & test_labels == 1)) - ...
        mean(Y_tests(test_p_inc_vec == p_inc_test & test_labels == 0)), ...
        test_p_inc_opts, 'uni', 1);

    a_v(:,cnt_save) = a; 
    b_v(:,cnt_save) = b; 
    t_v(:,cnt_save) = i; 

    bar_alphaW_v(:,cnt_save) = sum(reshaped_W(opts.selrow,:),'all') / sum(W(:)); 
    
    cnt_save = cnt_save + 1;
end

%% Save


res = struct; 
res.t = t_v;
res.a = a_v;
res.b = b_v;

res.bar_alpha_W = bar_alphaW_v - bar_alphaW_0; 
res.W = W_s;  

dY_objvsnoise_test = reshape(dY_objvsnoise_test, [1, size(dY_objvsnoise_test)]);
res.dY_objvsnoise_test = dY_objvsnoise_test;

res.vec_train_completeness = vec_train_completeness(t_v);
res.points_to_record_W = points_to_record_W; 

end

