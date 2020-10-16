function varargout = return_spikecount_matrix(spike_times, bin_width, bin_limit, bin_edges)

spike_times = check_spike_time_cell(spike_times); 

if ~exist('bin_edges', 'var') || any(isnan(bin_edges)) || isempty(bin_edges)
    if ~exist('bin_limit', 'var') || isnan(bin_limit) || isempty(bin_limit)
        low_lim = minnonnan(cellfun(@minnonnan, spike_times, 'uni', 1));
        upp_lim = maxnonnan(cellfun(@maxnonnan, spike_times, 'uni', 1));
        bin_limit = [low_lim, upp_lim];
    end
    
    if length(bin_limit) ~= 2 || bin_limit(1) >= bin_limit(2)
        error('The input of bin limits must be a vector of 2 elements of increasing order');
    end
    
    bin_limit(2) = bin_limit(1) + bin_width * ceil(diff(bin_limit)/bin_width);
    bin_edges = bin_limit(1):bin_width:bin_limit(2);
end

spike_counts = cellfun(@(x) ...
    to_row_vec(histcounts(x(:), 'BinEdges', bin_edges)), ...
    spike_times, 'uni', 0);
spike_counts = vertcat(spike_counts{:}); 

bin_centers = (bin_edges(1:end-1) + bin_edges(2:end))/2; 
results = {spike_counts, bin_centers, bin_edges}; 

varargout = results(1:nargout);


end

function x = minnonnan(x)
x = min(x, [], 'all', 'omitnan');
end

function x = maxnonnan(x)
x = max(x, [], 'all', 'omitnan');
end