function [slope, thres] = mahoncharpier(input, activity, slope, thres, eta_ip, ~)
dthres = eta_ip * tanh(activity - input); 
thres = thres + dthres; 
end