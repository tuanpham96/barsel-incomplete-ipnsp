function W = oja_bilinear(input, output, W, eta_sp_ltp, eta_sp_ltd)
dW = eta_sp_ltp * (output * input' - (output.^2) .* W) - eta_sp_ltd * (output + input');
W = W + dW; 
end
