function d = calculate_Hamming_distance(X, Y)
[X, Y] = check_same_length(X, Y); 
d = sum(xor(X, Y)); 
end

% benchmark for different methods
% L = 1e5; 
% N = 1e3; 
% 
% X = logical(randi(2, L, 1) - 1);
% Y = logical(randi(2, L, 1) - 1);
% 
% t0s = nan(N,1);
% t1s = nan(N,1);
% t2s = nan(N,1); 
% for i = 1:N
%     % slowest
%    ti = tic; abc = sum(X == true & Y == false) + sum(X == false & Y == true); 
%    t0s(i) = toc(ti);    
%    
%    % medium
%    ti = tic; abc = sum(X + Y == 1); 
%    t1s(i) = toc(ti); 
%    
%    % fatest
%    ti = tic; abc = sum(xor(X, Y)); 
%    t2s(i) = toc(ti);
% end
% figure; hold on; 
% histogram(t0s, 'displayname', 'sumsum'); 
% histogram(t1s, 'displayname', 'sumcheck'); 
% histogram(t2s, 'displayname', 'sumxor');