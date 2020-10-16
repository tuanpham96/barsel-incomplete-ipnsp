function poiss_result = stationary_poisson(X, poiss_param, x_param)

dt = poiss_param.dt; 
tau_ref = return_field_value(poiss_param, 'tau_ref', NaN); 

x_identity = upper(x_param.identity); 

bin_width = x_param.bin_width; 
bin_limit = return_field_value(x_param, 'bin_limit', NaN); 
bin_edges = return_field_value(x_param, 'bin_edges', NaN); 
return_type = return_field_value(x_param, 'return_type', 'count');

if contains(x_identity, {'RATE', 'COUNT'})
    if contains(x_identity, 'COUNT'), X = X/bin_width; end % convert to rate if x is count 
else 
    error('The identity of X needs to contain "rate" or "count", which "%s" does not', x_identity);
end 

r0 = to_col_vec(mean(X, 2)); 
[num_neurons, num_bins] = size(X); 
t = 0:dt:(num_bins*bin_width);
num_tstep = length(t); 

if isnan(tau_ref) 
    poiss_spkt = r0*dt > rand(num_neurons, num_tstep); 
else 
    poiss_spkt = zeros(num_neurons, num_tstep); 
    r = r0; 
    for i = 1:num_tstep
        p_spk = r*dt > rand(num_neurons, 1); 
        poiss_spkt(p_spk,i) = 1; 
        r(p_spk) = 0; 
        
        non_spk = ~p_spk; 
        dr = dt*(r0 - r)/tau_ref; 
        r(non_spk) = r(non_spk) + dr(non_spk);
    end
end 

spike_times = cell(num_neurons, 1); 
for i = 1:num_neurons 
    spike_times{i} = t(logical(poiss_spkt(i,:))); 
end

poiss_count = return_spikecount_matrix(spike_times, bin_width, bin_limit, bin_edges);

switch upper(return_type) 
    case 'COUNT', poiss_result = poiss_count; 
    case 'RATE', poiss_result = poiss_count/bin_width; 
    case {'TIMES', 'TIME'}, poiss_result = spike_times; 
end 

end
