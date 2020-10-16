function res = dsigmoid_dslope(inp, slope, thres) 
res = (inp-thres) .* exp(-slope.*(inp-thres)) ./ ...
    (exp(-slope.*(inp-thres)) + 1).^2;
end