function structobj = sample_struct(bin_ind, max_numbins, ...
    spike_times, spike_count, bin_edges)

obj_binrange = bin_ind:(bin_ind + max_numbins - 1);
obj_binedges = bin_edges(obj_binrange(1):obj_binrange(end)+1);
obj_timebounds = obj_binedges([1,end]);
obj_spikecount = spike_count(:, obj_binrange);
obj_spiketimes = sample_spiketimes(spike_times, obj_timebounds);

structobj = struct('count', obj_spikecount, 'time', {obj_spiketimes}, ...
    'bin_edges', obj_binedges);
        
end