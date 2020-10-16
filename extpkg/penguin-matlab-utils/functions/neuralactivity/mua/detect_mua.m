function [t_events, recordings, description] = detect_mua(recordings, time_info, detect_props) 
% recordings: filtered recordings "n_channels x n_times" 
% time_info: time in seconds or sampling frequency Fs 
% detect_props: struct with fields: 
%   threshold_scale: scaling for sigma_noise
%   min_event_distance: minimum event distance, in seconds 
%   detect_sign: (-1) to detect negatives, (+1) to detect positives,
%       (0,nan,[]) default to look for absolutes 
%   max_wave_duration: maximum wave duration around the, in seconds
%   padding_ratio: padding ratio of non-waves, default = 0 

% Reference on threshold: 
% Quiroga04: Quiroga et al, 2004, Neural Comp. Unsupervised Spike Detection and Sorting with Wavelets and Superparamagnetic Clustering

factor_abs_med = 0.6745; % refer to Quiroga04

[threshold_scale, min_event_distance_time, detect_sign, max_wave_duration_time, padding_ratio] = ...
    return_field_value(detect_props, ...
    'threshold_scale', 4, 'min_event_distance', 2e-3, 'detect_sign', 0, ...
    'max_wave_duration', nan, 'padding_ratio', 0); 

if length([threshold_scale,min_event_distance_time]) ~= 2, error('The threshold and distance must be scalar'); end 

if isempty(max_wave_duration_time), max_wave_duration_time = 3e-3; end % meaning extract wave, with default 3ms 
extract_wave = ~isnan(max_wave_duration_time); 
if extract_wave && nargout < 2, error('Wave extraction must be accompanied by 2 second output argument'); end 

abs_recordings = abs(recordings); 
if isvector(abs_recordings), abs_recordings = to_row_vec(abs_recordings); end
[num_channels, num_times] = size(abs_recordings);
if isscalar(time_info), time_info = (0:num_times-1)*time_info; end 

dtime = time_info(2) - time_info(1); 
min_event_distance_point = ceil(min_event_distance_time/dtime); 
if extract_wave, wave_window = ones(ceil(max_wave_duration_time/dtime),1); end 

t_events = cell(num_channels, 1);

sigma_noise = median(abs_recordings/factor_abs_med, 2);
mua_thres = threshold_scale * sigma_noise;

if detect_sign ~= 0 && ~isnan(detect_sign) && ~isempty(detect_sign)
    abs_recordings = sign(detect_sign) * recordings; 
end 

for chan_ind = 1:num_channels
    abs_rec_chan = abs_recordings(chan_ind,:);
    [~, loc_events_chan] = findpeaks(abs_rec_chan, ...
        'MinPeakHeight', mua_thres(chan_ind), 'MinPeakDistance', min_event_distance_point);
    t_events{chan_ind} = time_info(loc_events_chan); 
    if extract_wave
        event_vec = zeros(length(time_info), 1); 
        event_vec(loc_events_chan) = 1; 
        event_vec = conv(event_vec, wave_window, 'same') > 0;
        if padding_ratio == 0 % save some time when assign directly 
            recordings(chan_ind, ~event_vec) = 0;
        else 
            recordings(chan_ind, ~event_vec) = recordings(chan_ind, ~event_vec) * padding_ratio;
        end
    end
end

if num_channels == 1, t_events = t_events{1}; end 

description = sprintf('detect: (threshold_scale = %.2f) (min_event_distance = %.2f ms) ', threshold_scale, min_event_distance_time*1000); 
if extract_wave, description = sprintf('%s (max_wave_duration = %.2f ms)', description, max_wave_duration_time); end
end