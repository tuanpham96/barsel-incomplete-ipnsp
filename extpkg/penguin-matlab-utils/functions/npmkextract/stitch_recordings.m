function stitch_recordings(stitch_params) 

data_drive = return_field_value(stitch_params, 'data_drive', '.'); 

[file_condition, file_analysis, max_window_per_file] = return_field_value(...
    stitch_params, 'condition', '', 'analysis', '', 'max_window_per_file', [], '*strict*', true); 

switch upper(regexprep(file_analysis, {'_', '-'}, ''))
    case 'LFP'
        file_list = find_datanames(data_drive, 'condition', file_condition, ...
            'analysis', 'lfp', 'type', 'mat'); 
         stitch_lfp(file_list, data_drive, max_window_per_file); 
    case {'PRESPIKEMDA', 'PRESPIKESMDA', 'MOUNTAINSORT'} 
        file_info = find_datanames(data_drive, 'condition', file_condition, ...
            'analysis', 'prespikes', 'type', 'mat');
        file_list = find_datanames(data_drive, 'condition', file_condition, ...
            'analysis', 'prespikes', 'type', 'mda');
        [file_info, file_list] = matchfile_info_recordings(file_info, file_list);
        if length(file_info) ~= length(file_list) 
            error('The prespike stitching requires similar number of ''.mat'' and ''.mda'' files'); 
        end 
        stitch_prespikesmda(file_list, file_info, data_drive, max_window_per_file);
    case {'PRESPIKEDAT', 'PRESPIKESDAT', 'PRESPIKEBIN', 'PRESPIKESBIN', 'KILOSORT'} 
        error('Currently unimplemented yet, please contact if needed.');
    otherwise
        error('The file analysis type "%s" is not allowed', file_analysis); 
end 


end

function stitch_lfp(file_list, data_drive, max_window_per_file)
num_files = length(file_list); 
if num_files == 0, error('No files to stitch'); end

[DATA_INFO, num_stichedfile, stichsections] = ...
    obtain_datainfo(file_list, max_window_per_file);

load(file_list(1).path, 'dt'); 
load(file_list(1).path, 'Fs'); 

fullfile_prefix = file_list(1).identifier; 
fullfile_prefix(regexp(fullfile_prefix, 'lfp'):end) = '';
fullfile_prefix = fullfile(data_drive, fullfile_prefix); 

cur_tic = tic; 
fprintf('Stiching %s (total = %d) \n', fullfile_prefix, num_stichedfile); 
for sec_ind = 1:num_stichedfile 
    section = bound_minmax(stichsections(sec_ind,:), 1, num_files); 
    loaded_data = arrayfun(@(fn) load(fn.path, 'recordings', 'time'), ...
        file_list(section(1):section(end)), 'uni', 1);
    recordings = horzcat(loaded_data.recordings); 
    time = vertcat(loaded_data.time); 
    DATA_INFO.section = struct('num', sec_ind, 'range', time([1,end])); 
    saved_file_name = sprintf('%slfpstiched_%03d.mat', fullfile_prefix, sec_ind-1);
    save(saved_file_name, 'DATA_INFO', 'recordings', 'time', 'dt', 'Fs');
    
    fprintf('\t+ %s done. Current time elpased = %.2f sec.\n', saved_file_name, toc(cur_tic));
end
end

function stitch_prespikesmda(file_list, file_info, data_drive, max_window_per_file)
num_files = length(file_list); 
if num_files == 0, error('No files to stitch'); end

[DATA_INFO, num_stichedfile, stichsections] = ...
    obtain_datainfo(file_info, max_window_per_file);

load(file_info(1).path, 'dt'); 
load(file_info(1).path, 'Fs'); 

fullfile_prefix = file_list(1).identifier; 
fullfile_prefix(regexp(fullfile_prefix, 'prespikes'):end) = '';
fullfile_prefix = fullfile(data_drive, fullfile_prefix); 

cur_tic = tic; 
fprintf('Stiching %s (total = %d) \n', fullfile_prefix, num_stichedfile); 

for sec_ind = 1:num_stichedfile 
    section = bound_minmax(stichsections(sec_ind,:), 1, num_files); 
    fileind2read = section(1):section(end); 
    loaded_time = arrayfun(@(fn) load(fn.path, 'time'), file_info(fileind2read), 'uni', 1);
    time = vertcat(loaded_time.time); 
    DATA_INFO.section = struct('num', sec_ind, 'range', time([1,end])); 
    save(sprintf('%sprespikesstiched_%03d.mat', fullfile_prefix, sec_ind-1), ...
        'DATA_INFO', 'time', 'dt', 'Fs');
    
    loaded_mda = arrayfun(@(fn) readmda(fn.path), file_list(fileind2read), 'uni', 0); 
    recordings = horzcat(loaded_mda{:});
    saved_fullname = sprintf('%sprespikesstiched_%03d.mda', fullfile_prefix, sec_ind-1); 
    writemda(int16(recordings),saved_fullname, 'int16');
    
    fprintf('\t+ %s done. Current time elpased = %.2f sec.\n', saved_fullname, toc(cur_tic));
  
end

end

function [DATA_INFO, num_stichedfile, stichsections] = ...
    obtain_datainfo(file_list, max_window_per_file)
num_splitfiles = length(file_list); 
load(file_list(1).path, 'DATA_INFO'); 
split_window = DATA_INFO.split_window; 
num_splits_per_stichedfile = floor(max_window_per_file/split_window); 
[num_stichedfile, stichsections] = split_range(1, num_splitfiles, num_splits_per_stichedfile, 1);
DATA_INFO.section = struct('num', [], 'range', []); 
DATA_INFO.split_window = max_window_per_file;
end