function extract_accelorometer(DATA_INFO, file_suffix, select_accel_chan, select_qual_chan, split_window, smooth_window)

if ~exist('select_qual_chan', 'var')
    error('Must input also file suffix, selection for accelerometer and transmission quality channels');
end

if ~exist('split_window', 'var'), split_window = 'full'; end
if ~exist('smooth_window', 'var'), smooth_window = 1; end

accel_filepath = [DATA_INFO.fullfile_prefix, file_suffix];
accel_header = openNSx(accel_filepath, 'noread');
accel_header = accel_header.MetaTags;

total_num_sec = double(accel_header.DataPointsSec);
Fs = double(accel_header.SamplingFreq);
dt = 1/Fs;

[num_section, sections] = split_range(0, total_num_sec, split_window, dt);

cur_tic = tic; 
fprintf('Accelorometer extracting (total = %d) ... \n', num_section); 
for i = 1:num_section
    time_range = sections(i,:);
    identifier = sprintf('%03d',i-1);
    accel_struct = openNSx(accel_filepath, 'read', 'duration', time_range, 'sec');
    
    accel = double(accel_struct.Data(select_accel_chan,:));
    accel = sqrt(sum(accel.^2,1))';
    accel = smoothdata(accel, 'gauss', ceil(1/dt));
    
    trans_qual = double(accel_struct.Data(select_qual_chan,:))';
    trans_qual(trans_qual < 0) = 0;
    trans_qual(trans_qual > 0) = 1;
    
    time = (0:(length(accel)-1))' * dt + time_range(1);
    
    description = sprintf([...
        'ACCEL %s \n' ...
        'Preprocess: square -> sum -> sqrt -> smooth (Gauss)\n' ...
        'dt [s]: scalar \n' ...
        'Fs [Hz]: scalar \n' ...
        'accel [1]: [n_times x 1] \n' ...
        'trans_qual [bool]: [n_times x 1] \n' ...
        'time [s]: [n_times x 1]'], identifier);
    
    DATA_INFO.section = struct('num', i-1, 'range', time_range);
    DATA_INFO.split_window = split_window;
    DATA_INFO.smooth_window = smooth_window;
    DATA_INFO.select_accel_chan = select_accel_chan;
    DATA_INFO.select_qual_chan = select_qual_chan;
    
    save(sprintf('%s_accel_%s.mat', DATA_INFO.saved_prefix, identifier), ...
        'DATA_INFO', 'dt', 'Fs', 'description', ...
        'accel', 'trans_qual', 'time');
    
    fprintf('\t + accel-%s file extracted. Current elapsed time was %.2f seconds \n',identifier, toc(cur_tic));
end


fprintf('Accelorometer file extracted. Current elapsed time was %.2f seconds \n', toc(cur_tic)); 
end
