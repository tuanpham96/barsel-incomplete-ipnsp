function s = calculate_SMC(X,Y)
[X, Y] = check_same_length(X, Y); 
s = 1 - sum(xor(X, Y))/length(X); 
end

% L = 1e5; 
% N = 1e4; 
% 
% X = logical(randi(2, L, 1) - 1);
% Y = logical(randi(2, L, 1) - 1);
% 
% t0s = nan(N,1);
% t1s = nan(N,1);
% t2s = nan(N,1); 
% for i = 1:N
%     % slowest
%    ti = tic; abc = (sum(X == true & Y == true) + sum(X == false & Y == false))/length(X); 
%    t0s(i) = toc(ti);    
%    
%    % slightly faster than the next
%    ti = tic; abc = 1 - sum(xor(X, Y))/length(X); 
%    t1s(i) = toc(ti); 
%    
%    % fast-ish
%    ti = tic; abc = sum(xor(X, ~Y))/length(X); 
%    t2s(i) = toc(ti);
% end
% 
% figure; hold on; 
% histogram(t0s, 'displayname', 'sumsum', 'normalization', 'pdf'); 
% histogram(t1s, 'displayname', '1-sumxor', 'normalization', 'pdf'); 
% histogram(t2s, 'displayname', 'sumxornot', 'normalization', 'pdf');
% legend('show')