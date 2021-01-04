kr_v = 0:0.01:2;
a_v = arrayfun(@(kr) randfact2alpha(kr), kr_v, 'uni', 0); 
a_v = vertcat(a_v{:});

figure; hold on;
plot(kr_v,a_v)
xlabel('$k_r$');
ylabel('$\alpha$');
legend('$\alpha_{L1}$', '$\alpha_{L2}$');
function a = randfact2alpha(kr)
a = zeros(1000,2);
for i = 1:1000
    w = zeros(10);
    w(8,:) = 1;
    w = w(:)' + kr*rand(1,100);
    w = w / sqrt(sum(w.^2));
    w_r = reshape(w,[10,10]);
    a(i,:) = [sum(w_r(8,:)) / sum(w_r(:)), ...
        sqrt(sum(w_r(8,:).^2)) / sqrt(sum(w_r(:).^2))] ;
end
a = mean(a); 
end

