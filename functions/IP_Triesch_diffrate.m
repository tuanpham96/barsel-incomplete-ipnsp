function [a, b] = IP_Triesch_diffrate(X, Y, a, b, eta_ip_a, eta_ip_b, mu_Y) 
common_factor = 1 - (2+1./mu_Y).*Y + (Y.^2)./mu_Y;  
da = eta_ip_a * (1./a + X .* common_factor);
db = eta_ip_b * common_factor;
a = a + da; 
b = b + db; 
end