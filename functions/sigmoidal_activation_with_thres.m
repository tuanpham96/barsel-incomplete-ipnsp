function y = sigmoidal_activation_with_thres(x, a, theta, p)
b = log(1-p) - log(p); 
y = 1./(1 + exp(-a.*x + a.*theta + b)); 
end