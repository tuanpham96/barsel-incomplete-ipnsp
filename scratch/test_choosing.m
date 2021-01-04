% [desired_a0_mat, desired_b0_mat] = meshgrid(unique(params_tbl.a_init), unique(params_tbl.b_init));

max_ncolors = 5;

figure; hold on; 
cmap = return_colorbrewer('Set1',max_ncolors); 
colormap(cmap); 
mat = nan(5,5); 
image(mat', 'AlphaData', ~isnan(mat)');
caxis([1,size(cmap,1)]);
set(gcf, 'KeyPressFcn', @(~,~) choose_pixel(gcf));

function choose_pixel(fig)
    if ~strcmpi(get(fig, 'currentcharacter'), 'c')
        return; 
    end
    
    ax = findall(gcf, 'type', 'axes');
    mat = get(ax, 'children').CData';
    
    [x_coord,y_coord] = ginput(1); 
    x_coord = min(max(round(x_coord), 1), size(mat,1));
    y_coord = min(max(round(y_coord), 1), size(mat,2));
    
    mat_val = mat(x_coord, y_coord);
    
    max_mat = 0; 
    if any(~isnan(mat(:)))
        max_mat = max(mat(:), [], 'omitnan');
    end
    
    if isnan(mat_val) 
        mat_val = max_mat+1;
    else
        mat_val = nan;
    end
    
    mat(x_coord, y_coord) = mat_val;

    cla;  
    image(ax, mat', 'AlphaData', ~isnan(mat'));
    
 
end