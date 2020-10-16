function [slope, thres] = triesch_logistic(input, output, thres, slope, eta_ip, mu_output) 
[a, b] = triesch_og(input, output, slope, -slope.*thres, eta_ip, mu_output); 
slope = a; 
thres = -b./a; 
end