function v = normalize_minmax(v, low, upp)
switch nargin 
    case 1
        low = 0; 
        upp = 1;
    case 2
        error('Either not enter neither lower or upper bounds, or enter both. The default bound is [0, 1]');
end

min_v = min(v(:)); 
max_v = max(v(:)); 

v = (v - min_v) ./ (max_v - min_v); 
v = (upp - low) * v + low; 
end