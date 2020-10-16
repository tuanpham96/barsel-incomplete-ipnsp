function m = replacenan(m,v)
if nargin == 1, v = 0; end
m(isnan(m)) = v; 
end