function [waves, phases, envelopes, time] = waves_and_hilbertphase(recordings, Fs, wave_names) 
load('lfp_bands.mat', 'lfp_bands'); 

if ~exist('wave_names', 'var')    
    band_fields = fieldnames(lfp_bands);
    wave_names = band_fields(~contains(band_fields, 'ratio', 'IgnoreCase', true));
end

time = (0:(size(recordings,2)-1))*1/Fs; 
if iscell(wave_names) 
    phases = empty_struct_like(wave_names); 
    waves = empty_struct_like(wave_names); 
    envelopes = empty_struct_like(wave_names); 
    for i = 1:length(wave_names) 
        wave_name = wave_names{i}; 
        [w,p,e,~] = waves_and_hilbertphase(recordings, Fs, wave_name);
        waves.(wave_name) = w; 
        phases.(wave_name) = p;
        envelopes.(wave_name) = e; 
    end
    
    return;
end

freq_cutoffs = lfp_bands.(wave_names)/(Fs/2); 

if any(freq_cutoffs == 0)    
    [b, a] = butter(4, freq_cutoffs(2), 'low');
elseif any(freq_cutoffs == 1)
    [b, a] = butter(4, freq_cutoffs(1), 'high');
else
    [b, a] = butter(4, freq_cutoffs, 'bandpass');
end
waves = mean(filtfilt(b, a, recordings'),2);
hibert_wave = hilbert(waves); 
envelopes = abs(hibert_wave); 
phases = atan2(imag(hibert_wave),real(hibert_wave));
end