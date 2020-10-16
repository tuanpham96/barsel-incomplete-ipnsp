function res = discretize_colorlevels(v, range_v, color_var)
% discretize_colorlevels(v, range_v, num_colors) 
% discretize_colorlevels(v, range_v, cmap_matrix)
% INPUT
% v:            value vector 
% range_v:      2-element vector of hard [min, max] bound 
% color_var:    either a scalar to specify the number of levels 
%               or a "num_levels x 3" matrix of the colormap values 
% RETURN
% res:          if `color_var` is scalar -> return a list of integers to
%               specify the level 
%               if `color_var` is a matrix -> will return a "length(v) x 3"
%               matrix specify the colors for each value in v 

    if ~isvector(v), error('The first input has to be a vector'); end 
    if range_v(2) <= range(1), error('The second input needs to be a range of increasing values'); end 
    
    return_ind = false; 
    if length(color_var) == 1
        num_levels = color_var; 
        return_ind = true; 
    else
        if length(size(color_var)) ~= 2 || size(color_var, 2) ~= 3
            error('The third input can only be a scalar or N x 3 matrix'); 
        end 
        num_levels = size(color_var, 1);
    end 
    
    res = bound_minmax(ceil(num_levels*(v - range_v(1))/diff(range_v)),1,num_levels);
    if ~return_ind
        res = color_var(res, :); 
    end 
end