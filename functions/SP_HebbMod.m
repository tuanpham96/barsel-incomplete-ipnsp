function W = SP_HebbMod(X, Y, W, eta_sp_ltp, eta_sp_ltd)
dW = eta_sp_ltp * (Y * X') - eta_sp_ltd * (Y.^2 + (X').^2);
W = W + dW;
end