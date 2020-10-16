function s = calculate_Jaccard_similarity(X, Y)
[X, Y] = check_same_length(X, Y); 
s =  sum(X & Y)/sum(X | Y);
end


% L = 1e5;
% N = 1e3;
% 
% X = logical(randi(2, L, 1) - 1);
% Y = logical(randi(2, L, 1) - 1);
% 
% t0s = nan(N,1);
% t1s = nan(N,1);
% for i = 1:N
%     % slow
%     ti = tic; abc = sum(X == true & Y == true)/(length(X) - sum(X == false & Y == false));
%     t0s(i) = toc(ti);
%     
%     % fast
%     ti = tic; abc = sum(X & Y)/sum(X | Y);
%     t1s(i) = toc(ti);
%     
% end
% 
% figure; hold on;
% histogram(t0s, 'displayname', 'sumdivsum');
% histogram(t1s, 'displayname', 'sumlogical');
% legend('show')