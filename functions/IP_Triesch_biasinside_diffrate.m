function [a, b] = IP_Triesch_biasinside_diffrate(X, Y, a, b, eta_ip_a, eta_ip_b, mu_Y) 
common_factor = 1 - (2+1./mu_Y).*Y + (Y.^2)./mu_Y;  
da = eta_ip_a * (1./a + (X + b).* common_factor);
db = eta_ip_b * a .* common_factor;
a = a + da; 
b = b + db; 
end