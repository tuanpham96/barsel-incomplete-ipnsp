function W = oja_biquadratic(input, output, W, eta_sp_ltp, eta_sp_ltd)
dW = eta_sp_ltp * (output * input' - (output.^2) .* W) - eta_sp_ltd * (output.^2 + (input').^2);
W = W + dW; 
end