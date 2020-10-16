function res = calculate_LFP_powers_1channel(stft_chan, num_datapoints, freq_vec) 
load('lfp_bands.mat', 'lfp_bands'); 

band_fields = fieldnames(lfp_bands); 
ratio_fieldinds = contains(band_fields, 'ratio', 'IgnoreCase', true); 
ratio_fields = band_fields(ratio_fieldinds); 
nonratio_fields = band_fields(~ratio_fieldinds); 

num_ratio_elements = 2; % in the value 

pow_chan = stft_chan .*  conj(stft_chan) / num_datapoints; 
d_freq = freq_vec(2) - freq_vec(1); 

res = empty_struct_like(lfp_bands); 

for i = 1:length(nonratio_fields)
    band_name = nonratio_fields{i}; 
    freq_range = lfp_bands.(band_name); 
    check_frequency_range(freq_range); 
    res.(band_name) = d_freq * sum(pow_chan(freq_vec >= freq_range(1) & freq_vec <= freq_range(2),:), 1);  
end

for i = 1:length(ratio_fields) 
    ratio_name = ratio_fields{i};
    ratio_bands = lfp_bands.(ratio_name); 
    
    if length(ratio_bands) ~= num_ratio_elements
        error('There needs to be %d elements in the field %s', num_ratio_elements, ratio_bands);
    end
    
    band_pows = cell(num_ratio_elements, 1);
    
    for j = 1:num_ratio_elements
        band_element = ratio_bands{j};
        if ischar(band_element)
            if ~isfield(res, band_element), error('There is no band "%s".', band_element); end
            if isempty(res.(band_element)), error('The band "%s" has not been calculated.', band_element); end
            band_pows{j} = res.(band_element);
        else % needs to be a frequency range
            check_frequency_range(band_element);
            band_pows{j} = d_freq * sum(pow_chan(freq_vec >= band_element(1) & freq_vec <= band_element(2),:), 1);
        end
    end
    
    res.(ratio_name) = band_pows{1} ./ band_pows{2};
end

end


function check_frequency_range(freq_range)
error_mes = 'The frequency range needs to have 2 increasing numbers';  
if length(freq_range) ~= 2 || ~isnumeric(freq_range) 
    error(error_mes);
end 
if freq_range(1) >= freq_range(2) 
    error(erro_mes);
end
end