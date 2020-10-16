function visualize_colorbrewer(src) 
if ~exist('src', 'var'), src = 'plot'; end

switch upper(src) 
    case 'PNG'
        png_file = 'colorbrewer.png';
        if ~exist(png_file, 'file')
            error('There is no %s file', png_file);
        end
        [I,map] = imread(png_file);
        figure; imshow(imresize(I,0.2), map, 'InitialMagnification', 50);
        set(gca, 'position', [0,0,1,1]);
            
    case 'PLOT'

        load('colorbrewer.mat', 'colorbrewer');
        available_cmaps = fieldnames(colorbrewer);
        num_cmaps = length(available_cmaps);
        nrows = 10; ncols = ceil(num_cmaps/nrows);
        
        figure('color', 'w');
        for i = 1:num_cmaps
            cmap_name = available_cmaps{i};
            [cmap_mat, cmap_type] = return_colorbrewer(cmap_name);
            pseudo_values = linspace(-1, 1, size(cmap_mat,1));
            subplot(nrows, ncols, i);
            image(pseudo_values, 'cdatamapping', 'scaled');
            colormap(gca, cmap_mat);
            set(gca, 'xcolor', 'none', 'ycolor', 'none');
            title(sprintf('%s (%s)', cmap_name, cmap_type), 'fontsize', 12);
            
        end
    otherwise 
        error('%s is not allowed as an option of source', src);
end

end