function W = oja_og(input, output, W, eta_sp_ltp, ~)
dW = eta_sp_ltp * (output * input' - (output.^2) .* W);
W = W + dW; 
end
