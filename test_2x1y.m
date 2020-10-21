function res = test_2x1y(p_x, a, b, eta, T)
w1 = 1;
w2 = 1; 
w =  sqrt(w1^2 + w2^2);
w1 = w1 / w;
w2 = w2 / w;

alphaL1_v = zeros(T,1);
alphaL2_v = zeros(T,1);
w_v = zeros(T,2);
for t = 1:T

    x1 = rand < p_x; 
    x2 = 1 - x1;
    
    x = w1 * x1 + w2 * x2;  
    y = sigmoidal_activation(x, a, b); 
    
    w =  sqrt(w1^2 + w2^2);
    
    alphaL1_v(t) = w1/ (w1 + w2); 
    alphaL2_v(t) = w1/ w;
    
    w_v(t,:) = [w1,w2];
    
    dw1 = eta * y * x1; 
    dw2 = eta * y * x2; 
    w1 = w1 + dw1;
    w2 = w2 + dw2;
    
    
    w1 = w1 / w;
    w2 = w2 / w;

end

res.alpha_L1 = alphaL1_v;
res.alpha_L2 = alphaL2_v;
res.w1 = w_v(:,1);
res.w2 = w_v(:,2);
end