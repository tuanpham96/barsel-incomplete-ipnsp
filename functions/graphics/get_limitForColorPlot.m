function bound_v = get_limitForColorPlot(v, factor)
if ~exist('factor', 'var'), factor = 0.5; end 
v = sort(v(:)); 
bound_v = [v(1)-(v(2)-v(1))*factor, v(end)+(v(end)-v(end-1))*factor];
end