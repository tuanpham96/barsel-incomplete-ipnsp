function d = calculate_Jaccard_distance(X, Y)
d = 1 - calculate_Jaccard_similarity(X, Y); 
end

% L = 1e5;
% N = 2e4;
% 
% X = logical(randi(2, L, 1) - 1);
% Y = logical(randi(2, L, 1) - 1);
% 
% t0s = nan(N,1);
% t1s = nan(N,1);
% for i = 1:N
%     % slow
%     ti = tic; abc = sum(xor(X, Y))/sum(X | Y);
%     t0s(i) = toc(ti);
%     
%     % faster by a bit
%     ti = tic; abc = 1 - sum(X & Y)/sum(X | Y);
%     t1s(i) = toc(ti);
%     
% end
% 
% figure; hold on;
% histogram(t0s, 'displayname', 'sumxor', 'Normalization', 'cdf');
% histogram(t1s, 'displayname', '1-sim',  'Normalization', 'cdf');
% legend('show')