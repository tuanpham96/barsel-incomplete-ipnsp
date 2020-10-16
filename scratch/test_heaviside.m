T = 2e4; 
y_v = zeros(T,1);
b_v = zeros(T,1); 

mu_y = 0.1; 
eta_ip = 1e-3; 
inp_mu = 1; 
inp_sigma = 1; 

a = 1e-10; 
delta_fun = @(x) exp(-(x/a).^2) / (abs(a)*sqrt(pi));
d_delta = @(x) delta_fun(x) * (-2*x) / (pi * a^4); 

b = 0; 

inp_fun = @(T) inp_sigma * (randn(T, 1) + inp_mu); 

x_v = inp_fun(T); 

for t = 1:T
    x = x_v(t); 
    y = x + b >= 0; 
    
    delta_xb = delta_fun(x + b); 
    d_delta_xb = d_delta(x + b); 
    db = eta_ip * (d_delta_xb / max(delta_xb,a) - delta_xb/mu_y);
    b = b + db; 
    
    y_v(t) = y;
    b_v(t) = b;
end

figure; 
subplot(221); hold on; 
plot(y_v); 

subplot(222); hold on; 
plot(y_v);
% 
% subplot(223); hold on; 
% y_eval = y_v(T/2:T); 
% histogram(y_eval,100);
% % title(sprintf('$\\mu = %g, E(y) = %.4f$', mu_y, mean(y_eval))); 






