clc; clear; % close all;
run start_up.m; 

N_side = 10; 
N_X = N_side^2;
N_Y = 1;

a_Y = 5; 
b_Y = -1; 

% a_Y = 6; 
% b_Y = -5; 

mu_Y = 1/(2*N_side); 
eta_ip = 1e-2;
eta_sp = 1e-3; 

num_train = 20e4;

p_bar = 1/(2*N_side); 

% W_YX = normalize_weight(lognrnd(-1,0.5,N_Y, N_X));
% W_YX = normalize_weight(rand(N_Y, N_X));
% W_YX = normalize_weight(betarnd(2,8,[N_Y, N_X]));
% W_YX = normalize_weight(randn(N_Y, N_X));

W_YX = zeros(N_side); 
W_YX(8,:) = 1; 
W_YX = W_YX(:); 
W_YX = normalize_weight(W_YX' + rand(N_Y,N_X));

% X_trains = arrayfun(@(~) generate_bar_input(N_side, p_bar), 1:num_train, 'uni', 0); 

num_train_complete = ceil(9/10 * num_train);
num_train_incomplete = num_train - num_train_complete; 
X_trains = [arrayfun(@(~) generate_bar_input(N_side, p_bar),1:num_train_complete, 'uni', 0), ...
            arrayfun(@(~) generate_bar_input_incomplete(N_side, p_bar,0.1), 1:num_train_incomplete, 'uni', 0)];

Y_trains = zeros(N_Y,num_train); 

W_summary = zeros(N_X,num_train); 
a_summary = zeros(1,num_train); 
b_summary = zeros(1,num_train); 
pre_Ys = zeros(1,num_train); 
alphaW_summary = zeros(1, num_train);
tic
for i = 1:num_train
    
    % train
    X_train = X_trains{i}; 
  
    pre_Y_train = W_YX * X_train ; 
    Y_train = sigmoidal_activation(pre_Y_train, a_Y, b_Y);
  
    Y_trains(:,i) = Y_train;
    pre_Ys(i) = pre_Y_train;
    W_summary(:,i) = W_YX;
    a_summary(i) = a_Y;
    b_summary(i) = b_Y;
    
    w_reshape = normalize_minmax(reshape(W_YX, [N_side,N_side])); 
    alphaW_summary(i) = sum(w_reshape(8,:))/sum(w_reshape(:));
    % update
    [a_Y, b_Y] = IP_Triesch(pre_Y_train, Y_train, a_Y, b_Y, eta_ip, mu_Y);
    W_YX = normalize_weight(SP_Hebbian(X_train, Y_train, W_YX, eta_sp));

end
toc

%% Plot example inputs 
num_inp_toplot = 25; 
Us_toplot = X_trains(randperm(num_train, num_inp_toplot)); 
nrows = ceil(sqrt(num_inp_toplot)); ncols = nrows; 
figure('position', [0,0,0.5,1]); colormap('gray'); 
for i = 1:num_inp_toplot
    subplot(nrows,ncols,i); hold on;
    image(reshape(Us_toplot{i}, [N_side,N_side])); 
    
    if i == 1
        title({'example inputs', sprintf('$p_{bar}=%g$', p_bar)}); 
    end
    
    set(gca, 'xtick', '', 'ytick', '', 'box', 'on'); 
    axis square;
end
%%

n_splt = 16; 
nrows = ceil(sqrt(n_splt)); ncols = nrows; 
n_factor_plt = round(num_train/n_splt); 

figure; colormap('gray');


for i = 1:n_splt
    subplot(nrows,ncols,i); hold on; 
    image_with_strict_limits(reshape(W_summary(:,(i-1)*n_factor_plt + 1), [N_side,N_side]));

    axis square; 
end
%%

Ys = zeros(N_side, 2); 
for i = 1:N_side
    X = zeros(N_side); 
    X(i,:) = 1; X = X(:);
    Ys(i,1) = sigmoidal_activation(W_YX * X, a_Y, b_Y);
    
    X = zeros(N_side); 
    X(:,i) = 1; X = X(:);
    Ys(i,2) = sigmoidal_activation(W_YX * X, a_Y, b_Y);
end

figure; 
subplot(231); hold on; plot(Ys);

subplot(232); hold on; plot(Y_trains, 'k.'); 
subplot(233); hold on;
histogram(W_summary(:,1))
histogram(W_summary(:,num_train * 0.3))
histogram(W_summary(:,num_train * 0.6))
histogram(W_summary(:,end))
legend('init', '1/3', '2/3', 'end'); 

% subplot(234); hold on; histogram(Y_trains(:,end/2:end)); 
% title(sprintf('$\\mu=%g, r=%.5f$', mu_Y, mean(Y_trains(:,end/2:end))));

subplot(235); 
hold on; plot(a_summary); 
yyaxis right; plot(b_summary);

subplot(236); hold on; 
plot(W_summary(randperm(N_X, 10),1:100:end)');