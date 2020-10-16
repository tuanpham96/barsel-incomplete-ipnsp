function W = SP_Cov(X, Y, eY, W, eta_sp)
dW = eta_sp * ((Y - eY)* X');
W = W + dW; 
end
