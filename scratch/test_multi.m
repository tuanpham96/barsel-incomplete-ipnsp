clc; clear; close all;
run start_up.m;

save_progress = false; 

N = 20;
N_U = N^2;
N_Y = 100;

% a = randn(N_Y,1);
% b = randn(N_Y,1);
a = 5*rand(N_Y,1); 
b = 5*(2*randn(N_Y,1) - 1);
% a = 2; 
% b = -3;
mu = 1/(6*N);
% mu = 0.1;

eta_ip = 1e-2;
eta_sp = 1e-2;

num_train = 20e4;

p_bar = 1/(6*N);
W_init = normalize_weight(rand(N_Y, N_U) -0.5);
% W_init = normalize_weight(rand(N_Y, N_U));

W = W_init;

% U_trains = arrayfun(@(~) generate_bar_input(N, p_bar), 1:num_train, 'uni', 0);

tic
U_trains = arrayfun(@(~) generate_dir_input(N, p_bar), 1:num_train, 'uni', 0);
toc

tic
for i = 1:num_train
    
    U = U_trains{i};
    
    X = W * U ;
    Y = sigmoidal_activation(X, a, b);
    
    % update
    [a, b] = IP_Triesch(X, Y, a, b, eta_ip, mu);
%     W = normalize_weight(SP_Hebbian(U, Y, W, eta_sp));
    W = normalize_weight(SP_HebbMod(U, Y, W, eta_sp));

end

toc

W_final = W;

%%
Ws_reshape_size = [10,10]; 
pad_size = 2;
Ws = arrayfun(@(i) padarray(reshape(W_final(i,:), [N,N]), pad_size * [1,1], nan, 'post'), 1:N_Y, 'uni', 0); 
Ws = cell2mat(reshape(Ws, Ws_reshape_size));

cmap = return_colorbrewer('RdBu_R', 50); 

figure; colormap(cmap); 
image(Ws, 'AlphaData', ~isnan(Ws));
caxis([-1,1]*max(abs(Ws(:))));
axis square; colorbar;
set(gca, 'visible', 'off');

% figure; 
% colormap('gray');
% nrows = 10; ncols = 10;
% for i = 1:N_Y
%     subplot(nrows, ncols, i); hold on; 
%     image_with_strict_limits(reshape(W_final(i,:), [N,N])); 
%     axis square;
%     set(gca, 'visible', 'off'); 
% end

%%
% N_Y = 100; 
% tic
% Ws = arrayfun(@(~) test_W, 1:N_Y, 'uni', 0);
% toc
% 
% 
% figure; 
% colormap('gray');
% nrows = 10; ncols = 10;
% for i = 1:N_Y
%     subplot(nrows, ncols, i); hold on; 
%     image_with_strict_limits(Ws{i});
%     axis square;
%     set(gca, 'visible', 'off'); 
% end

%%
% figure; 
% N_U_2plt = 25; 
% for i = 1:N_U_2plt
%     subplot(5,5,i); hold on;
%     Ui = reshape(U_trains{i}, [N,N]); 
%     image_with_strict_limits(Ui);
%     axis square;
%     set(gca, 'visible', 'off'); 
% end
