function W = SP_Hebbian(X, Y, W, eta_sp)
dW = eta_sp * (Y * X');
W = W + dW; 
end
