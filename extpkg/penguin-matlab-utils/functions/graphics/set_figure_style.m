function set_figure_style(varargin)
first_inp = varargin{1}; 
curr_argind = 2; 
if strcmp(first_inp.Type, 'figure')
    fig = first_inp; 
else
    fig = gcf; 
    curr_argind = 1; 
end 
normz_pos = varargin{curr_argind}; curr_argind = curr_argind + 1;   
scale     = varargin{curr_argind}; 

set_figure_style_with_param(fig, normz_pos, scale); 
end

function set_figure_style_with_param(fig, normz_pos, scale)
set(fig, 'unit', 'normalized', 'position', normz_pos) ;
set(fig, ...
    'PaperUnits', 'centimeters', ...
    'PaperSize', scale*[1,normz_pos(4)/normz_pos(3)],...
    'PaperPosition', [0,0,1,normz_pos(4)/normz_pos(3)]*scale,...
    'InvertHardcopy', 'off', ...
    'Color', 'w');
end