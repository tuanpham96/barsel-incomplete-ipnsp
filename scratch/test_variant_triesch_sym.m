syms a theta x y mu
y = 1 / (1 + exp(-a*(x-theta))); 
%%
dy_dx = diff(y,x);
dy_da = diff(y,a); 
dy_dtheta = diff(y,theta); 
%%
simplify(dy_da - (x-theta) * y * (1-y))
simplify(dy_dtheta - (-a) * y * (1-y))
simplify(dy_dx - a * y * (1-y))
%%
D = -log(dy_dx) + y/mu;
c = 1 - (2+1/mu)*y + y^2 / mu; 

simplify(diff(D, a) - (-1/a -(x-theta)*c))
simplify(diff(D, theta) - a*c)


%%
syms a b theta x y mu p
b = log(1-p) - log(p); 
y = 1 / (1 + exp(-a*(x-theta) + b)); 
D = -log(diff(y, x)) + y/mu;
c = 1 - (2+1/mu)*y + y^2 / mu; 

simplify(diff(D, a) - (-1/a -(x-theta)*c))
simplify(diff(D, theta) - a*c)