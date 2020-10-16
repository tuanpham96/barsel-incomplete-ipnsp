function extract_recordings(DATA_INFO, file_suffix, resample_struct, select_record_chan, split_window)
% TODO: restructure, when have time, to have less repeated codes 

%% Parse input 
if ~exist('resample_struct', 'var')
    error('Must input also file suffix and resampling struct');
end

[lfp_resample, preSpikes_resample, mua_resample] = return_field_value(resample_struct, ...
    'lfp', [], 'pre_spikes', [], 'mua', []); 
if isempty(lfp_resample) && isempty(preSpikes_resample) && isempty(mua_resample)
    error('The resampling struct needs to have either "lfp" or "pre_spikes" to continue extraction'); 
end 

% Parse LFP extraction/resampling struct 
extract_LFP = false; 
if ~isempty(lfp_resample) 
    extract_LFP = true; 
    [lfp_Fs, lfp_filterband, lfp_filtertype, lfp_filteroder] = return_field_value(lfp_resample, ...
        'Fs', nan, 'band', nan, 'type', nan, 'order', 3); 
    if xor(all(isnan(lfp_filterband)), isnan(lfp_filtertype))
        error('The LFP resample struct needs to have both or neither filter band and type'); 
    end
end

% Parse preprocessed recordings for spikes extraction/resampling struct 
extract_preSpikes = false; 
if ~isempty(preSpikes_resample) 
    extract_preSpikes = true; 
    [prespikes_Fs, prespikes_filterband, prespikes_filtertype, prespikes_filterorder, prespikes_saveformat] = ...
        return_field_value(preSpikes_resample, ...
        'Fs', nan, 'band', nan, 'type', nan, 'order', 3, 'saveformat', 'dat'); 
    if xor(all(isnan(prespikes_filterband)), isnan(prespikes_filtertype))
        error('The spike resample struct needs to have both or neither filter band and type'); 
    end
end

% Parse MUA extraction/resampling struct 
extract_MUA = false; 
extract_MUAwaves = false; 
if ~isempty(mua_resample) 
    extract_MUA = true; 
    [mua_Fs, mua_filterband, mua_filtertype, mua_filterorder, mua_detectprops, mua_wavesaveformat] = ...
        return_field_value(mua_resample, ...
        'Fs', nan, 'band', nan, 'type', nan, 'order', 4, 'detect', struct, 'saveformat', 'dat'); 
    if xor(all(isnan(mua_filterband)), isnan(mua_filtertype))
        error('The MUA resample struct needs to have both or neither filter band and type'); 
    end
    if isempty(mua_detectprops) 
        error('The MUA detection property struct must not be empty for double-checking purposes'); 
    end
    if isfield(mua_detectprops, 'max_wave_duration'), extract_MUAwaves = true; end 
end

% Information specific to data 
if ~exist('select_record_chan', 'var'), select_record_chan = 'full'; end 
if ~exist('split_window', 'var'), split_window = 'full'; end

%% Read header and split 
rec_filepath = [DATA_INFO.fullfile_prefix, file_suffix];
rec_header = openNSx(rec_filepath, 'noread');
rec_header = rec_header.MetaTags;

total_num_sec = double(rec_header.DataPointsSec);
Fs_original = double(rec_header.SamplingFreq);

[num_section, sections] = split_range(0, total_num_sec, split_window);
DATA_INFO.Fs_original = Fs_original; 

%% Extract 
cur_tic = tic; 
fprintf('Recordings extracting (total = %d) ... \n', num_section); 
for section_ind = 1:num_section
    Fs = Fs_original; dt = 1/Fs; 
    time_range = sections(section_ind,:);
    identifier = sprintf('%03d',section_ind-1);
    original_recordings = openNSx(rec_filepath, 'read', ...
        'channels', select_record_chan, 'duration', time_range, 'sec');
    original_recordings = original_recordings.Data;
    
    %% LFP extraction 
    if extract_LFP
        recordings = double(original_recordings);
        
        resampled_text = '';
        if ~isnan(lfp_Fs)
            recordings = resample(recordings', lfp_Fs, Fs_original)';
            Fs = lfp_Fs;
            dt = 1/Fs;
            resampled_text = sprintf('resampled (%d to %d)', Fs_original, lfp_Fs);
        end
        
        filter_text = '';
        if ~isnan(lfp_filterband)
            if strcmpi(lfp_filtertype, 'notch')
                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',lfp_filterband-1,'HalfPowerFrequency2',lfp_filterband+1, ...
                    'DesignMethod','butter','SampleRate',Fs);
                recordings = filtfilt(d, recordings')';
            else
                [b, a] = butter(lfp_filteroder, lfp_filterband/(Fs/2), lfp_filtertype);
                recordings = filtfilt(b, a, recordings')';
            end
            filter_text = sprintf('-> %s filter %s Hz ', lfp_filtertype, sprintf('%.1f ', lfp_filterband));
        end
        
        time = (0:(size(recordings, 2)-1))' * dt + time_range(1);
        
        description = sprintf([...
            'Recordings LFP %s \n' ...
            'Preprocess: %s %s \n' ...
            'dt [s]: scalar \n' ...
            'Fs [Hz]: scalar \n' ...
            'recordings [1]: [n_channels x n_times(_resampled)] \n' ...
            'time [s]: [n_times x 1]'], identifier, resampled_text, filter_text);
        
        DATA_INFO.preprocess = lfp_resample;
        DATA_INFO.section = struct('num', section_ind-1, 'range', time_range);
        DATA_INFO.split_window = split_window;
        DATA_INFO.select_record_chan = select_record_chan;
        
        save(sprintf('%s_lfp_%s.mat', DATA_INFO.saved_prefix, identifier), ...
            'DATA_INFO', 'dt', 'Fs', 'description', ...
            'recordings', 'time');
        clear time recordings 
        fprintf('\t + lfp-%s file extracted. Current elapsed time was %.2f seconds \n',identifier, toc(cur_tic));
    end
    
    %% preSpikes 
    if extract_preSpikes
        recordings = double(original_recordings);
        
        resampled_text = '';
        if ~isnan(prespikes_Fs)
            recordings = resample(recordings', prespikes_Fs, Fs_original)';
            Fs = prespikes_Fs;
            dt = 1/Fs;
            resampled_text = sprintf('resampled (%d Hz to %d Hz)', Fs_original, prespikes_Fs);
        end
        
        filter_text = '';
        if ~isnan(prespikes_filterband)
            if strcmpi(prespikes_filtertype, 'notch')
                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',prespikes_filterband-1,'HalfPowerFrequency2',prespikes_filterband+1, ...
                    'DesignMethod','butter','SampleRate',Fs);
                recordings = filtfilt(d, recordings')';
            else
                [b, a] = butter(prespikes_filterorder, prespikes_filterband/(Fs/2), prespikes_filtertype);
                recordings = filtfilt(b, a, recordings')';
            end
            [b, a] = butter(prespikes_filterorder, prespikes_filterband/(Fs/2), prespikes_filtertype);
            recordings = filtfilt(b, a, recordings')';
            filter_text = sprintf('-> %s filter %s Hz', prespikes_filtertype, sprintf('%.1f ', prespikes_filterband));
        end
        
        time = (0:(size(recordings, 2)-1))' * dt + time_range(1);
        
        description = sprintf([...
            'Recordings preprocessed for spikes extraction %s \n' ...
            'Preprocess: %s %s \n' ...
            'dt [s]: scalar \n' ...
            'Fs [Hz]: scalar \n' ...
            'recordings [1]: [n_channels x n_times(_resampled)] \n' ...
            'time [s]: [n_times x 1]'], identifier, resampled_text, filter_text);
        
        
        DATA_INFO.preprocess = preSpikes_resample;
        DATA_INFO.section = struct('num', section_ind-1, 'range', time_range);
        DATA_INFO.split_window = split_window;
        DATA_INFO.select_record_chan = select_record_chan;
        
        % MATLAB file to save general information 
        saved_fileprefix = sprintf('%s_prespikes_%s', DATA_INFO.saved_prefix, identifier); 
        save([saved_fileprefix '.mat'], 'DATA_INFO', 'dt', 'Fs', 'description', 'time');
        
        % Binary file to save the recording matrix
        saved_fullname = [saved_fileprefix '.' lower(prespikes_saveformat)];
        switch upper(prespikes_saveformat)
            case {'DAT', 'BIN'}
                fid = fopen(saved_fullname, 'w');
                fwrite(fid, int16(recordings), 'int16');
                fclose(fid);
            case 'MDA'
                writemda(int16(recordings),saved_fullname, 'int16');
            otherwise
                error('Only DAT, BIN or MDA is acceptable for saved binary format');
        end
        clear time recordings
        fprintf('\t + prespikes-%s file extracted. Current elapsed time was %.2f seconds \n',identifier, toc(cur_tic));
        
    end
    
    %% MUA events 
    if extract_MUA
        recordings = double(original_recordings);
        
        resampled_text = '';
        if ~isnan(mua_Fs)
            recordings = resample(recordings', mua_Fs, Fs_original)';
            Fs = mua_Fs;
            dt = 1/Fs;
            resampled_text = sprintf('resampled (%d to %d)', Fs_original, mua_Fs);
        end
        
        filter_text = '';
        if ~isnan(mua_filterband)
            if strcmpi(mua_filtertype, 'notch')
                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',mua_filterband-1,'HalfPowerFrequency2',mua_filterband+1, ...
                    'DesignMethod','butter','SampleRate',Fs);
                recordings = filtfilt(d, recordings')';
            else
                [b, a] = butter(mua_filterorder, mua_filterband/(Fs/2), mua_filtertype);
                recordings = filtfilt(b, a, recordings')';
            end
            filter_text = sprintf('-> %s filter %s Hz ', mua_filtertype, sprintf('%.1f ', mua_filterband));
        end
        time = (0:(size(recordings, 2)-1))' * dt + time_range(1);
        
        [event_times, recordings, detect_text] = detect_mua(recordings, time, mua_detectprops); 
        detect_text = sprintf('-> %s', detect_text);
        
        description = sprintf([...
            'MUA events extraction %s \n' ...
            'Preprocess: %s %s %s\n' ...
            'dt [s]: scalar \n' ...
            'Fs [Hz]: scalar \n' ...
            '(recordings [1]: [n_channels x n_times(_resampled)]) \n' ...
            'event_times [s]: {n_channels x 1}[n_events_in_channel x 1] \n' ...
            'time [s]: [n_times x 1]'], identifier, resampled_text, filter_text, detect_text);
        
        DATA_INFO.preprocess = mua_resample;
        DATA_INFO.section = struct('num', section_ind-1, 'range', time_range);
        DATA_INFO.split_window = split_window;
        DATA_INFO.select_record_chan = select_record_chan;
        
        % MATLAB file to save general information and MUA event times         
        saved_fileprefix = sprintf('%s_mua_%s', DATA_INFO.saved_prefix, identifier); 
        save([saved_fileprefix '.mat'], 'DATA_INFO', 'dt', 'Fs', ...
            'description', 'event_times', 'time');         
        
        % Binary file to save the recording matrix 
        if extract_MUAwaves
            saved_fullname = [saved_fileprefix '.' lower(mua_wavesaveformat)]; 
            switch upper(mua_wavesaveformat) 
                case {'DAT', 'BIN'} 
                    fid = fopen(saved_fullname, 'w');
                    fwrite(fid, int16(recordings), 'int16');
                    fclose(fid);
                case 'MDA'                    
                    writemda(int16(recordings),saved_fullname, 'int16');
                otherwise
                    error('Only DAT, BIN or MDA is acceptable for saved binary format');
            end
        end
                
        clear time recordings            
        fprintf('\t + mua-%s file extracted. Current elapsed time was %.2f seconds \n',identifier, toc(cur_tic));
    end
    
end

fprintf('Recordings file(s) extracted. Current elapsed time was %.2f seconds \n', toc(cur_tic)); 
end
