function spike_times = check_spike_time_cell(spike_times)

if ~iscell(spike_times)
    spike_times = {spike_times};
end

if ~all(cellfun(@isnumeric, spike_times, 'uni', 1))
    error('The `spike_times` needs to be a cell array, each element of which is a vector of numeric numbers');
end

end