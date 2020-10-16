clc; clear; close all;
run start_up.m; 

% rng(1022);
N_side = 10; 
N_X = N_side^2;
N_Y = 1;

a_Y = 5; 
b_Y = -10; 

mu_Y = 1/(2*N_side); 
mu_Y = 0.05;
eta_ip = 1e-2; 
eta_sp = 1e-3; 

num_train = 20e4;

p_bar = 1/(2*N_side); 
% W_YX = normalize_weight(lognrnd(-1,0.5,N_Y, N_X));
W_YX = normalize_weight(rand(N_Y, N_X));

% W_YX = normalize_weight(betarnd(2,8,[N_Y, N_X]));
% W_YX = normalize_weight(randn(N_Y, N_X));

X_trains = arrayfun(@(~) generate_bar_input(N_side, p_bar), 1:num_train, 'uni', 0); 
Y_trains = zeros(N_Y,num_train); 

W_summary = zeros(N_X,num_train); 
a_summary = zeros(1,num_train); 
b_summary = zeros(1,num_train); 
pre_Ys = zeros(1,num_train); 

num_avr = 10; 
cnt_avr = 0; 

Da = 0; 
Db = 0;
DW = 0; 

tic
for i = 1:num_train
    
% X = X / sum(X); 
    % train
    X_train = X_trains{i}; 
    
% norm_X = sqrt(sum(X_train.^2));
% if norm_X > 0, norm_X = 1/norm_X;end
% X_train  = X_train * norm_X;
    pre_Y_train = W_YX * X_train ;
    Y_train = 1 ./ (1 + exp( - (a_Y * pre_Y_train + b_Y))); 
%     Y_train = sigmoidal_activation(pre_Y_train, a_Y, b_Y);
  
    Y_trains(:,i) = Y_train;
    pre_Ys(i) = pre_Y_train;
    W_summary(:,i) = W_YX;
    a_summary(i) = a_Y;
    b_summary(i) = b_Y;
    
    % update
%     [a_Y, b_Y] = IP_Triesch(pre_Y_train, Y_train, a_Y, b_Y, eta_ip, mu_Y);
%     W_YX = (SP_Hebbian(X_train, Y_train, W_YX, eta_sp));

common_factor = 1 - (2+1./mu_Y).*Y_train + (Y_train.^2)./mu_Y;  
da = eta_ip * (1./a_Y +  pre_Y_train .* common_factor);
db = eta_ip * common_factor;

a_Y = a_Y + da; 
b_Y = b_Y + db; 


% dW = eta_sp * (Y_train * X_train') - Y_train.^2 .* W_YX;
% W_YX = (W_YX + dW); 

% dW = eta_sp * (Y_train * X_train') - eta_sp/20 * (Y_train.^2 + X_train'.^2);
% W_YX = normalize_weight(W_YX + dW); 
 
% dW = eta_sp * (Y_train * X_train');
% W_YX = normalize_weight(W_YX + dW); 

dW = eta_sp * (Y_train * X_train');
W_YX = normalize_weight(W_YX + dW); 

% if cnt_avr == num_avr
% a_Y = a_Y + Da/num_avr; 
% b_Y = b_Y + Db/num_avr; 
% % W_YX = normalize_weight(W_YX + DW/num_avr); 
% cnt_avr = 0;
% Da = 0;
% Db = 0; 
% % DW = 0; 
% else
%     Da = Da + da; 
%     Db = Db + db; 
% %     DW = DW + dW; 
%     cnt_avr = cnt_avr + 1; 
% end
end
toc

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