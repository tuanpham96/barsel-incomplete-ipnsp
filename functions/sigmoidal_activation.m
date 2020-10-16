function y = sigmoidal_activation(x, a, b)
y = 1./(1 + exp(-(a.*x + b))); 
end