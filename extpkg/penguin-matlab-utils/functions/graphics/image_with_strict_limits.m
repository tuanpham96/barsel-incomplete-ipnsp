function image_with_strict_limits(varargin) 
% image_with_strict_limits(mat)
% image_with_strict_limits(mat, clim)
% image_with_strict_limits(ax, mat)
% image_with_strict_limits(ax, mat, clim)

switch(nargin)
    case 1
        mat = varargin{1}; 
    case 2 
        if ismatrix(varargin{1}) && isnumeric(varargin{1})
            [mat, clim] = varargin{:};           
        elseif ishandle(varargin{1})
            [ax, mat] = varargin{:}; 
        else 
            error('Input types not allowed!');
        end
    case 3
        [ax, mat, clim] = varargin{:};
    otherwise 
        error('Number of input must not exceed 3');
end
if ~exist('ax', 'var'), ax = gca; end
if ~ismatrix(mat), error('Matrix is missing'); end
if ~ishandle(ax), error('Axes handle missing'); end 
if exist('clim', 'var'), caxis(ax, clim); end 

sz = size(mat); 
image(ax, mat', 'AlphaData', ~isnan(mat')); 
xlim(ax, [0.5, sz(1)+0.5]);
ylim(ax, [0.5, sz(2)+0.5]);

end