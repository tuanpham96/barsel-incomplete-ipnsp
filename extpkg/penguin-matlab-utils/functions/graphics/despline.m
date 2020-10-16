function despline(varargin)
% Usage:
%   despline(ax, offsetfactor)
%           ax:             axes handle 
%           offsetfactor:   can be a number or 2-element vector ([x,y])  
%   despline(ax);                   % offsetfactor = 1 
%   despline(offsetfactor)          % ax = gca 
%   despline([ax1, ax2, ...])
%   despline([ax1, ax2, ...], offsetfactor)
%   despline(fig_handle)
%   despline(fig_handle, offsetfactor)
%   despline('all')                 % 'all' is equil to all axes in gcf  
%   despline('all', offsetfactor)
% 
% Inpiration: offsetAxes by Anne Urai, 2016 
% https://www.mathworks.com/matlabcentral/fileexchange/52351-offsetaxes-ax
% 
% Modified by Tuan Pham, 2019 
% - Oct 06, 2019: fixes to adapt to log scale, add offsetfactor
% - Oct 25, 2019: fixes minor bugs 
% Note: 
% - For stable usage, set the limits and scales of the axes before using.
%   Minor ticks still unstably controlled (unnecessary minor ticks may
%   disappear after the first limit setting - meaning may not redrawing 
%   minor ticks smaller than lower bound of limit) 
% - Currently not including Z axis 
%% Parse input 
narginchk(0, 2); 

if nargin == 1 
    arg = varargin{1}; 
    if isnumeric(arg), offsetfactor = arg;
    else, ax = arg; 
    end
end

if nargin == 2
    ax = varargin{1}; 
    offsetfactor = varargin{2};
end

if ~exist('ax', 'var'), ax = gca; end
if ~exist('offsetfactor', 'var'), offsetfactor = 1; end 

if ischar(ax)
    if any(strcmpi(ax, {'all', 'allaxes', 'all-axes' 'all_axes'})), ax = gcf; 
    else, error('A character input could only be "ALL"'); end
end 

if isa(ax, 'matlab.ui.Figure')
    ax = findall(ax, 'type', 'axes');
end

if isvector(ax) && length(ax) > 1
    arrayfun(@(x) despline(x, offsetfactor), ax); 
    return; 
end

%% Execute 
if ~isa(ax, 'matlab.graphics.axis.Axes'), error('The first argument needs to be an axes handle'); end 
if ~isnumeric(offsetfactor), error('Second argument needs to be numeric'); end
offsetfactor_unit = 1/30; 
offsetfactor = offsetfactor * offsetfactor_unit; 
if length(offsetfactor) == 1, offsetfactor = [1,1]*offsetfactor; end
expandLowBound(ax, 'x', offsetfactor(1));
expandLowBound(ax, 'y', offsetfactor(2));

addlistener(ax, 'MarkedClean', @(obj,event) resetVertex(ax, offsetfactor));

end

%% Functions for resetting vertex data 
function resetVertex(ax, offsetfactor)

if any(strcmp({ax.XScale, ax.YScale}, 'log'))
    resetVertexLogScaleGeneral(ax, offsetfactor); 
else
    resetVertexLinearScale(ax); 
end

end


function resetVertexLinearScale(ax) 
if ~isempty(ax.XTick)
    ax.XRuler.Axle.VertexData(1,:) = get_bound(ax.XTick);
end
if ~isempty(ax.YTick)
    ax.YRuler.Axle.VertexData(2,:) = get_bound(ax.YTick);
end
end 

function resetVertexLogScaleGeneral(ax, offsetfactor) 
resetVertexLogScaleSpecific(ax, 'x', offsetfactor);
resetVertexLogScaleSpecific(ax, 'y', offsetfactor);
end 

function resetVertexLogScaleSpecific(ax, ax_type, offsetfactor)
ax_type = upper(ax_type); 
switch ax_type
    case 'X', ax_row = 1; 
    case 'Y', ax_row = 2; 
    otherwise, error('%s not accepted as an axis type, only "x" or "y"', ax_type);
end

ax_limit = ax.([ax_type 'Lim']); 
ax_tickbound = ax.([ax_type 'Tick']);

if isempty(ax_tickbound) 
    return;
else
    ax_tickbound = get_bound(ax_tickbound);
end

low_limit = ax_limit(1);
if strcmp(ax.([ax_type 'Scale']), 'linear')
    ax.([ax_type 'Ruler']).Axle.VertexData(ax_row,:) = (ax_tickbound - low_limit)/diff(ax_limit);
else 
    if low_limit <= 0 % special case: re-do low limit based on data
        ax.([ax_type 'Lim'])(1) = min(arrayfun(@(x) min(x.([ax_type 'Data'])), allchild(ax)));
        expandLowBound(ax, ax_type, offsetfactor(ax_row)); 
        ax_limit = ax.([ax_type 'Lim']); low_limit = ax_limit(1);
    end 
    % log scale: need to scale in log 
%     ax_tickbound(1) = max(ax_limit(1), ax_tickbound(1)); 
%     ax_tickbound(2) = max(ax_limit(2), ax_tickbound(2));
    ax.([ax_type 'Ruler']).Axle.VertexData(ax_row,:) = log10(ax_tickbound/low_limit)/diff(log10(ax_limit));
end

end

%% Functions for expanding low bounds
function expandLowBound(ax, ax_type, factor)
ax_type = upper(ax_type); 
if strcmp(ax.([ax_type 'Scale']), 'log')     
    ax.([ax_type 'Lim']) = logOffsetLowBound(ax.([ax_type 'Lim']), factor); 
else
    ax.([ax_type 'Lim']) = linearOffsetLowBound(ax.([ax_type 'Lim']), factor); 
end
ridUnecessaryMinorTicks(ax, ax_type); 
end

function ridUnecessaryMinorTicks(ax, ax_type)
ax_type = upper(ax_type); 
ax_ruler = ax.([ax_type 'Ruler']); 
ax_tickbound = get_bound(get(ax, [ax_type 'Tick'])); 
ax_ruler.MinorTickValues(ax_ruler.MinorTickValues < ax_tickbound(1) | ...
   ax_ruler.MinorTickValues > ax_tickbound(2)) = [];

end 
function new_limit = linearOffsetLowBound(axis_limit, factor)
offset_by = (axis_limit(2) - axis_limit(1))*factor; 
new_limit = axis_limit + [-1, 0]*offset_by;
end

function new_limit = logOffsetLowBound(axis_limit, factor)
if axis_limit(1) <= 0, error('Log scale cannot have non-positive lower limit'); end 
offset_by = log10(axis_limit(2)/axis_limit(1))*factor; 
new_limit = axis_limit .* 10.^([-1, 0]*offset_by);
end


