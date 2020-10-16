function mkwave2spiketimes(waveform_path, saved_prefix, minimum_ratings, timeunit_mkwave, padsec_betweenfiles) 
max_nums_displayed_on_line = 20; 

if isempty(waveform_path), waveform_path = uigetdir; end
if ~exist('saved_prefix', 'var'), saved_prefix = 'out'; end 
if ~exist('timeunit_mkwave', 'var'), timeunit_mkwave = 'ms'; end
if ~exist('padsec_betweenfiles', 'var'), padsec_betweenfiles = 5; end 

if waveform_path == 0, return; end
% This is because apparently in MK waveforms, the spiketimes time unit is
% in milisecond. This is just to convert everything into SECONDS 
switch upper(timeunit_mkwave) 
    case {'S', 'SEC', 'SECOND', 'SECONDS'}, timeunit_scale = 1; 
    case {'MS', 'MSEC', 'MILISECOND', 'MILISECONDS'}, timeunit_scale = 1/1000;
    otherwise, error('Invalid time units, only allowed seconds or miliseconds');
end

% Obtain files 
waveform_files = dir(fullfile(waveform_path, 'waveforms*.mat')); 
parent_folder = waveform_files(1).folder; 
[file_inds, ~] = listdlg('ListString', {waveform_files.name}, 'PromptString', ...
    'Select the files to be processed', 'ListSize', [300, 300]); 
waveform_files = waveform_files(file_inds); 
waveform_files = arrayfun(@(x) fullfile(parent_folder, x.name), waveform_files, 'uni', 0);

% Read information about source files
load(waveform_files{1}, 'waveforms'); 
sourceFiles = waveforms.sourceFiles.filenames;
num_sections = size(waveforms.sourceFiles.spikesByFile, 1);
    
% Start reading 
num_files = length(waveform_files); 
select_spiketimes = cell(num_files, 1); 
select_info = cell(num_files, 1); 
section_spiketimebounds = nan(num_files, num_sections); 

fprintf('Reading ... '); 
for file_ind = 1:num_files
    file_loc = waveform_files{file_ind}; 
    load(file_loc, 'waveforms'); 
    
    % Extract data 
    electrodeID = waveforms.electrode;    
    units = to_col_vec(waveforms.units); 
    spike_times = timeunit_scale * to_col_vec(waveforms.spikeTimes); 
    ratings = waveforms.ratings.ratings; 
    
    % Construct section by source files 
    section_spikeindices = waveforms.sourceFiles.spikesByFile;
    section_spiketimebounds(file_ind,:) = spike_times(section_spikeindices(:,2)); 
    
    if num_sections ~= size(section_spikeindices, 1), error('Mismatch number of source files'); end 
    sections = arrayfun(@(x) x*ones(diff(section_spikeindices(x,:))+1,1), 1:num_sections, 'uni', 0); 
    sections = vertcat(sections{:}); 
    
    % Selection and save to variables
    select_units = to_col_vec(find(ratings >= minimum_ratings));
    %     select_spiketimes{file_ind} = arrayfun(@(y) ...
    %         arrayfun(@(x) spiketimes(units == x & sections == y), ...
    %         select_units, 'uni', 0), 1:num_sections, 'uni', 0);
    
    select_spiketimes{file_ind} = arrayfun(@(y) ...
        arrayfun(@(x) spike_times(units == y & sections == x), ...
        1:num_sections, 'uni', 0), select_units, 'uni', 0);
    
    select_info{file_ind} = struct('file', file_loc, ...
        'electrode', electrodeID, 'units', num2cell(select_units)); 
    
    if mod(file_ind, max_nums_displayed_on_line) == 0, add_newline = '\n\t'; else, add_newline = ''; end
    fprintf([add_newline, '%d\t'], electrodeID); 
end
fprintf('\n ... Done reading. Start combining and saving ... \n'); 

select_spiketimes = vertcat(select_spiketimes{:}); 
info = vertcat(select_info{:}); 

if ~any(contains(saved_prefix, {'/','\'})), saved_prefix = fullfile(parent_folder, saved_prefix); end
    
    
if padsec_betweenfiles >= 0 && ~isnan(padsec_betweenfiles)
    % Stitching
    padded_maxdursec_between_sections = cumsum(max(section_spiketimebounds, [], 1) + padsec_betweenfiles);
    padding_between_sections = [0, padded_maxdursec_between_sections(1:end-1)];
    boundaries_between_sections = [padding_between_sections',padded_maxdursec_between_sections'];
    padding_between_sections = num2cell(padding_between_sections);
    spike_times = cellfun(@(x) cellfun(@(z,y) z+y, x, padding_between_sections, 'uni', 0),  ...
        select_spiketimes, 'uni', 0);
    spike_times = cellfun(@(x) vertcat(x{:}), spike_times, 'uni', 0);
    
    saved_file = [saved_prefix '_spiketimes.mat'];
    save(saved_file, 'spike_times', 'info', 'sourceFiles', ...
        'padsec_betweenfiles', 'boundaries_between_sections');
    fprintf(' Done! Saved as \n %s \n', saved_file);
else
    for section_ind = 1:num_sections 
        spike_times = cellfun(@(x) x{section_ind}, select_spiketimes, 'uni', 0); 
        saved_file = sprintf('%s_%03d_spiketimes.mat', saved_prefix, section_ind-1);
        sourceFile = sourceFiles{section_ind}; 
        save(saved_file, 'spike_times', 'info', 'sourceFile'); 
        fprintf(' + %s \n', saved_file);
    end
    fprintf('All done!\n'); 
end
end