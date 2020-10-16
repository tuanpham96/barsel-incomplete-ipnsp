function sampled = sample_spiketimes(spike_times, ref_bounds, offset)
low_bound = ref_bounds(1); 
upp_bound = ref_bounds(2); 

if ~exist('offset', 'var') 
    offset = 0;
elseif ischar(offset) 
    switch upper(offset)
        case {'LOW', 'LOWER', 'LOWER_BOUND'}
            offset = low_bound; 
        case {'UP', 'UPP', 'UPPER', 'UPPER_BOUND'}
            offset = upp_bound;
        otherwise 
            error('Only "low" or "up" are accepted');
    end
end

sampled = cellfun(@(x) x(x >= low_bound & x <= upp_bound) - offset, spike_times, 'uni', 0); 

end 