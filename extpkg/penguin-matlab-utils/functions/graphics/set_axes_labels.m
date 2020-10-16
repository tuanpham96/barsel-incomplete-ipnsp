function set_axes_labels(varargin)
first_inp = varargin{1}; 
curr_argind = 2;
if isa(first_inp, 'axes')
    ax = first_inp; 
else
    ax = gca; 
    curr_argind = 1; 
end 

xlbl = varargin{curr_argind}; curr_argind = curr_argind + 1;   
ylbl = varargin{curr_argind}; curr_argind = curr_argind + 1;   
ttl  = varargin{curr_argind};

set_axes_labels_with_param(ax, xlbl, ylbl, ttl); 
end

function set_axes_labels_with_param(ax, xlbl, ylbl, ttl) 
xlabel(ax, xlbl);
ylabel(ax, ylbl); 
title(ax, ttl); 
end