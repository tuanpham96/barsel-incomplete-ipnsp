function out = logistic_activation(inp, maxvalue, slope, thres) 
out = maxvalue./(1+exp(-slope.*(inp-thres))); 
end