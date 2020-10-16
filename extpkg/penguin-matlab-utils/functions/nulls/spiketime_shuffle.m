function shuf_res = spiketime_shuffle(spike_times, shuf_method, varargin)

spike_times = check_spike_time_cell(spike_times);

num_neurons = length(spike_times);

return_type = 'TIMES'; return_binwidth = nan; return_binedges = nan; 
if isstruct(shuf_method) 
    return_type = upper(return_field_value(shuf_method, 'return_type', 'times'));
    return_binwidth = return_field_value(shuf_method, 'return_bin_width', nan);
    return_binlimit = return_field_value(shuf_method, 'return_bin_limit', nan);   
    return_binedges = return_field_value(shuf_method, 'return_bin_edges', nan);   
    shuf_method = shuf_method.shuffle;
end

if ~contains(return_type, {'RATE', 'COUNT', 'TIME'}) 
    error('The return type "%s" is not allowed', return_type);
end

if contains(return_type, {'RATE', 'COUNT'}) && isnan(return_binwidth)
    error(['The return type "%s" needs to have an extra field "return_bin_width",' ...
        ' and an optional field "return_bin_limit" in "shuf_method'], return_type);
end

switch upper(shuf_method)
    case 'ISI', shuf_fun = @isi_shuffle;
    case 'DISPLACEMENT', shuf_fun = @displacement_shuffle;
    case 'EXCHANGE', shuf_fun = @spike_exchange_shuffle;
    otherwise, error('The method "%s" is not accepted', shuf_method);
end

spike_times = shuf_fun(spike_times, num_neurons, varargin{:});

if contains(return_type, 'TIME') 
    shuf_res = spike_times; 
    return;
end

shuf_res = return_spikecount_matrix(spike_times, return_binwidth, return_binlimit, return_binedges);

if contains(return_type, {'COUNT', 'RATE'})     
    if contains(return_type, 'RATE')
        shuf_res = shuf_res/return_binwidth;
    end 
end

end

function spike_times = isi_shuffle(spike_times, num_neurons, p_nrn)
if ~exist('p_nrn', 'var'), p_nrn = 1; end
for i = 1:num_neurons
    spkt = spike_times{i};
    num_spk = length(spkt);
    
    if rand > p_nrn || num_spk < 3, continue; end
    
    isi_nrn = diff(spkt);
    isi_nrn = isi_nrn(randperm(num_spk - 1));
    spike_times{i}(2:end) = spkt(1) + cumsum(isi_nrn);
end
end

function spike_times = displacement_shuffle(spike_times, num_neurons, max_dt, sign_on)
if ~exist('max_dt', 'var'), max_dt = 50e-3; end
if ~exist('sign_on', 'var'), sign_on = true; end

if sign_on
    low_dt = -max_dt; upp_minus_low_dt = 2*max_dt;
else
    low_dt = 0; upp_minus_low_dt = max_dt;
end

for i = 1:num_neurons
    spkt = spike_times{i};
    
    if isempty(spkt), continue; end 
    
    tmprl_displ = low_dt + upp_minus_low_dt*rand(size(spkt));
    new_spkt = spkt + tmprl_displ;
    % delete the spike times out of range of the original
    % for binning limit consistency -> maxtrix size consistency
    new_spkt(new_spkt < min(spkt) | new_spkt > max(spkt)) = [];
    spike_times{i} = new_spkt;
end
end

function spike_times = spike_exchange_shuffle(spike_times, num_neurons, p_pair_exch, p_spk_exch)
if ~exist('p_pair_exch', 'var'), p_pair_exch = 0.5; end
if ~exist('p_spk_exch', 'var'), p_spk_exch = 0.5; end

num_pair_exch = ceil(p_pair_exch*num_neurons*(num_neurons-1)/2);
for k = 1:num_pair_exch
    pair_ids = randperm(num_neurons, 2);
    spkt_i = spike_times{pair_ids(1)};
    spkt_j = spike_times{pair_ids(2)};
    
    num_spk_exch = ceil(p_spk_exch*min([length(spkt_i), length(spkt_j)]));
    
    if num_spk_exch < 1, continue; end
    
    exch_ind_i = randperm(length(spkt_i), num_spk_exch);
    exch_ind_j = randperm(length(spkt_j), num_spk_exch);
    
    exch_spk_i = spkt_i(exch_ind_i);
    exch_spk_j = spkt_j(exch_ind_j);
    
    spkt_i(exch_ind_i) = exch_spk_j;
    spike_times{pair_ids(1)} = spkt_i;
    
    spkt_j(exch_ind_j) = exch_spk_i;
    spike_times{pair_ids(2)} = spkt_j;
end
end