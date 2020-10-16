function [slope, thres] = mahoncharpier_activitygated(input, activity, slope, thres, eta_ip, ~)
dthres = eta_ip * tanh(activity .* (activity - input));  
thres = thres + dthres; 
end 