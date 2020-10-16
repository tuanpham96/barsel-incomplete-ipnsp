function extract_NEV2MDA(DATA_INFO, select_record_chan)

if ~exist('select_record_chan', 'var'), select_record_chan = 'full'; end 

fullfile_prefix = fullfile(DATA_INFO.data_drive, DATA_INFO.file_prefix);

fprintf('NEV file extracting ... \n'); 
nev = openNEV([fullfile_prefix '.nev'], 'nosave', 'nomat');
spikes_nev = nev.Data.Spikes; 
Fs = double(nev.MetaTags.SampleRes); 

dt = 1/Fs; 
nev_times = double(spikes_nev.TimeStamp)*dt;
electrode = double(spikes_nev.Electrode); 

if strcmpi(select_record_chan, 'full') 
    select_record_chan = 1:max(electrode);
end

spike_times = arrayfun(@(x) to_col_vec(double(nev_times(electrode == x))), select_record_chan', 'uni', 0);
waveforms = arrayfun(@(x) double(spikes_nev.Waveform(:, electrode == x)), select_record_chan', 'uni', 0);
time = (1:size(spikes_nev.Waveform,1))'*dt; 

description = sprintf([...
    'FROM NEV FILE \n' ...
    'dt [s]: scalar \n' ...
    'Fs [Hz]: scalar \n' ...
    'spike_times [s]: {n_channels x 1} -> [n_spikes x 1] \n' ...
    'waveforms [1]: {n_channels x 1} -> [n_time_points x n_waveforms] \n' ...
    'time [s]: [n_time_points x 1]']); 


DATA_INFO.select_record_chan = select_record_chan; 
save([DATA_INFO.saved_prefix '_nev.mat'], 'dt', 'Fs', 'description', 'DATA_INFO', ...
    'spike_times', 'waveforms', 'time'); 

fprintf('NEV file extracted. Elapsed time was %.2f seconds \n', toc); 
end

function recording = stich_waveforms_1channel(waveforms, spike_times, time_vec)
dt = time_vec(2) - time_vec(1); 
recording = zeros(length(time_vec), 1); 
if iscell(waveforms) 
    waveforms = vertcat(waveforms{:})'; 
end 
spkt_inds = round(spike_times/dt); 
[~, peak_inds] = max(waveforms, [], 1); 
wave_inds = repmat(1:size(waveforms,1), size(waveforms,2)); 
wave_inds = wave_inds - peak_inds + spkt_inds; 
recording(wave_inds(:)) = waveforms(:); 
end