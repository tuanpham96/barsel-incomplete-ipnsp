function y = sigmoidal_activation_biasinside(x, a, b)
y = 1./(1 + exp(-(a.*(x + b)))); 
end