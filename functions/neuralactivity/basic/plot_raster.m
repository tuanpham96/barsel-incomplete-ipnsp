function plot_raster(spike_times, raster_params)

if exist('raster_params', 'var')
    param_names = fieldnames(raster_params);
    for i = 1:length(param_names)
        prm = param_names{i};
        raster_params.(lower(prm)) = raster_params.(prm);
    end
else
    raster_params = struct;
end

[colors, markertype, linewidth, linelength, markersize] = ...
    return_field_value(raster_params, ...
    'color', [0,0,0], 'markertype', 'line', ...
    'linewidth', 1, 'linelength', 0.9, 'markersize', 10);

num_colors = size(colors, 1);
if num_colors > 1 && num_colors ~= length(spike_times)
    error('If use multiple colors, the number of colors have to match with length of `spike_times`');
end


if strcmpi(markertype, 'line')
    
    for i = 1:length(spike_times)
        spkt = spike_times{i};
        if isempty(spkt), continue; end
        clr = valid_color(colors, i);
        for j = 1:length(spkt)
            line([spkt(j),spkt(j)],[i,i+linelength], 'LineWidth', linewidth, 'Color', clr);
        end
        
    end
else
    hold on 
    for i = 1:length(spike_times)
        spkt = spike_times{i};
        if isempty(spkt), continue; end
        clr = valid_color(colors, i);
        plot(spkt, i*ones(length(spkt),1), 'LineStyle', 'none', 'Marker', markertype, ...
            'MarkerSize', markersize, 'MarkerEdgeColor', clr, ...
            'MarkerFaceColor', clr, 'Color', clr);
    end
    
end

end

function c = valid_color(colors, ind)
if size(colors, 1) > 1
    c = colors(ind, :);
else
    c = colors;
end
end