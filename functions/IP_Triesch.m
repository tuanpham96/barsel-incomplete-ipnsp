function [a, b] = IP_Triesch(X, Y, a, b, eta_ip, mu_Y) 
common_factor = 1 - (2+1./mu_Y).*Y + (Y.^2)./mu_Y;  
da = eta_ip * (1./a + X .* common_factor);
db = eta_ip * common_factor;
a = a + da; 
b = b + db; 
end
