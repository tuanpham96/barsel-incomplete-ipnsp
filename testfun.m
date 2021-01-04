function r = testfun(alpha,a,b,eta,x)
w = 1/sqrt(1 + ((1-alpha)/alpha)^2); 
y = sigmoidal_activation( w*x + ((1-alpha)/alpha)*w*(1-x), a, b); 
r = (w + eta * y .* (x - y*w)) ./ (w + eta * y .* (alpha - y*w));
end