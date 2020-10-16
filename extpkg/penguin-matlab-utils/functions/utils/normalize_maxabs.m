function v = normalize_maxabs(v, low, upp)
switch nargin 
    case 1
        low = 0; 
        upp = 1;
    case 2
        error('Either not enter neither lower or upper bounds, or enter both. The default bound is [0, 1]');
end

v = v / max(abs(v)); 
v = (upp - low) * v + low; 
end