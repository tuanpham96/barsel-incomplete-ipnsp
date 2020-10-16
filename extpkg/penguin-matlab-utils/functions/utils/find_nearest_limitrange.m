function [x_nearest, ind] = find_nearest_limitrange(v, x, max_look_range, return_first)
if ~exist('return_first', 'var'), return_first = 'val'; end

low = x - max_look_range; upp = x + max_look_range; 
rng2look = (v >= low & v <= upp);
first_ind_rng = find(rng2look, 1); 

[~, ind] = min(abs(v(rng2look)-x));
x_nearest = v(ind + first_ind_rng - 1);

if strcmpi(return_first, 'val') % do nothing 
elseif strcmpi(return_first, 'ind')
    [x_nearest, ind] = swap(x_nearest, ind); 
else 
    error('The parameter "return_first" can only be "val" or "ind" to indicate what to return first');
end
end