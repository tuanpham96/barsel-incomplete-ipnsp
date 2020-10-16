function hide_only_axis(varargin)
% Usage:
%   hide_only_axis          % meaning x,y,z will be hidden for gca
%   hide_only_axis('x')     % ax = gca 
%   hide_only_axis('xy')
%   hide_only_axis(ax)      % meaning x,y,z will be hidden for ax
%   hide_only_axis(ax, 'xyz')

narginchk(0, 2);
if nargin == 0
    ax = gca; 
    ax_types = 'xyz'; 
end

if nargin == 1
    if ischar(varargin{1})
        ax = gca; 
        ax_types = varargin{1};
    elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
        ax = varargin{1};
        ax_types = 'xyz';
    else
        error('Only an axes handle or a charater(array) representing x/y/z is accepted');
    end
end

if nargin == 2
    ax = varargin{1}; 
    ax_types = varargin{2};
end

arrayfun(@(c) hide_axis(ax, c), ax_types); 
end
function hide_axis(ax, ax_type)
% inpired from jonas
% https://www.mathworks.com/matlabcentral/answers/411646-how-to-make-axis-line-invisible-but-keep-xlabel-and-ylabel
ax_type = upper(ax_type);
ax_axes = ax.([ax_type 'Axis']);
for ind_ax = 1:length(ax_axes)
    ax_axis = ax_axes(ind_ax); 
    ax_color = ax_axis.Label.Color;
    ax_vison = strcmpi(ax_axis.Label.Visible, 'on');
    set(ax, [ax_type, 'Color'], 'none', [ax_type 'Tick'], '');
    
    if ax_vison
        ax_axis.Label.Color = ax_color;
        ax_axis.Label.Visible = 'on';
    end
end
end