function res = dsigmoid_dinp(inp, slope, thres) 
res = slope .* exp(-slope.*(inp-thres)) ./ ...
    (exp(-slope.*(inp-thres)) + 1).^2;
end
