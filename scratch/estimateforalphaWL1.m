hahaha = @(r,k,N) (1+r*0.5)./(1+N*0.5*r./k);

N = 100;
k_v = 1:99;
r_v = logspace(-2,2,100);
n = 100;

[r_v, k_v] = meshgrid(r_v,k_v);
alphas_Wres = zeros(size(r_v)); 

% for i = 1:size(r_v,1)
%     for j = 1:size(r_v,2)
%         alphas_Wres(i,j) = calc_alpha(r_v(i,j),k_v(i,j),N,n);
%     end
% end
%     
% theo_mean = hahaha(r_v,k_v,N);

calc_alpha(1,10,100,2e4)
calc_alpha(0.1,10,100,2e4)
calc_alpha(2,10,100,2e4)
calc_alpha(10,10,100,2e4)

function res = calc_alpha(r,k,N,n)
alpha_W = zeros(n,1);
for i = 1:n
    a = zeros(N,1);
    a(1:k)=1; 
    a = a + r*rand(size(a));
    a = a / sqrt(sum(a.^2));
    alpha_W(i) = sum(a(1:k))/sum(a(:)); 
end

res = mean(alpha_W);
end