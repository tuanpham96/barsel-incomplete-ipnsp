function out = sigmoid(inp, slope, thres) 
out = 1./(1+exp(-slope.*(inp-thres))); 
end