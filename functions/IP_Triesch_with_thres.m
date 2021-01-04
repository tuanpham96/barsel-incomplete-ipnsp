function [a, theta] = IP_Triesch_with_thres(X, Y, a, theta, eta_ip_a, eta_ip_theta, mu_Y) 
C       = 1 - (2+1./mu_Y).*Y + (Y.^2)./mu_Y; 

da      = 1./a + (X - theta) .* C;
dtheta  = -a .* C;

a       = a     + eta_ip_a      * da; 
theta   = theta + eta_ip_theta  * dtheta; 
end