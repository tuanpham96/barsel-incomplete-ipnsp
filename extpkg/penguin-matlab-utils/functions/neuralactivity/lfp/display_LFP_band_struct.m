function varargout = display_LFP_band_struct(band_name) 
load('lfp_bands.mat', 'lfp_bands') 
all_bands = fieldnames(lfp_bands); 
ratio_bands = all_bands(contains(all_bands, 'ratio')); 
actual_bands = all_bands(~contains(all_bands, 'ratio')); 

if nargin == 0 
    if nargout > 0, error('No output is available when displaying all'); end
else
    range_value = lfp_bands.(band_name); 
    if ~contains(band_name, 'ratio')
        band_range = get_range(range_value);
        band_name = get_label(band_name);
        varargout{1} = sprintf('%s: %s', band_name, band_range);
    else
        numer = get_range(range_value{1});
        denom = get_range(range_value{2});
        band_name = get_label(band_name);
        varargout{1} = sprintf('%s: [%s] / [%s]', band_name, numer, denom);
    end
    return; 
end
        

fprintf('The LFP bands are: \n'); 
for band_ind = 1:length(actual_bands) 
    band_name = actual_bands{band_ind}; 
    
    band_range = get_range(lfp_bands.(band_name));
    band_label = get_label(band_name); 
    fprintf('\t%s: %s \t (%s) \n', band_label, band_range, band_name); 
end
fprintf('The ratio measures are from: \n'); 
for ratio_ind = 1:length(ratio_bands) 
    ratio_name = ratio_bands{ratio_ind}; 
    ratio_vals = lfp_bands.(ratio_name); 
    numer = get_range(ratio_vals{1});
    denom = get_range(ratio_vals{2});
    ratio_label = get_label(ratio_name); 
    fprintf('\t%s: [%s] / [%s] \t (%s) \n', ratio_label, numer, denom, ratio_name); 
end
end
function x_range = get_range(x_range)
if isnumeric(x_range) 
    x_range = [vec2str(x_range, ' - ', ''), ' Hz'];
    x_range = regexprep(x_range, '.0 ', ' ');
elseif ischar(x_range)
    x_range = get_label(x_range);
end
end
function x_label = get_label(x_label)
x_label = strsplit(x_label, '_');
band_source = x_label{end};
x_label = strjoin(x_label(1:end-1), ' ');
x_label = sprintf('%s (%s)', x_label, band_source);
end