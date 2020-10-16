function [lfp_pows_avr,freq_vec,time_vec,lfp_pows_sep] = calculate_LFP_powers(recordings, Fs, T_window, T_overlap)
if ~exist('lfp_bands', 'var'), load('lfp_bands.mat', 'lfp_bands'); end

N_window = ceil(T_window*Fs); 
N_overlap = ceil(T_overlap*Fs); 
window_stft = hamming(N_window); 
num_channels = size(recordings, 1); 
[~, freq_vec, time_vec] = spectrogram(recordings(1,:), window_stft, N_overlap, [], Fs); 
freq_vec = to_col_vec(freq_vec); 
time_vec = to_col_vec(time_vec); 

clear lfp_pows_sep
lfp_pows_sep(num_channels, 1) = empty_struct_like(lfp_bands);
    
for i = 1:num_channels
    stft_chan = spectrogram(recordings(i,:), window_stft, N_overlap, [], Fs);     
    lfp_pows_sep(i) = calculate_LFP_powers_1channel(stft_chan, N_window, freq_vec); 
end

lfp_pows_sep = structarray_to_struct(lfp_pows_sep, 0); 
lfp_pows_sep = structfun(@(x) vertcat(x{:}), lfp_pows_sep, 'uni', 0);
lfp_pows_avr = structfun(@(x) mean(x, 1), lfp_pows_sep, 'uni', 0);


end

