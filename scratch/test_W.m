function W_final = test_W

N = 10;
N_U = N^2;
N_Y = 1;

% a = randn(N_Y,1);
% b = randn(N_Y,1);
a = 2;
b = -2;

mu = 1/(2*N);
% mu = 0.01; 

eta_ip = 1e-2;
eta_sp = 1e-2;

num_train = 10e4;

p_bar = 1/(2*N);
W_init = normalize_weight(rand(N_Y, N_U));

W = W_init;

U_trains = arrayfun(@(~) generate_bar_input(N, p_bar), 1:num_train, 'uni', 0);


for i = 1:num_train
    
    U = U_trains{i};
    
    X = W * U ;
    Y = sigmoidal_activation(X, a, b);
    
    % update
    [a, b] = IP_Triesch(X, Y, a, b, eta_ip, mu);
    W = normalize_weight(SP_Hebbian(U, Y, W, eta_sp));
    
end

W_final = reshape(W, [N,N]);
end