function res = read_neuron_file(file_names)
persistent recent_dir; 
if isempty(recent_dir), recent_dir = pwd; end 
if ~exist('file_names', 'var')
    sel_files = 'Select a set of files'; 
    sel_folder = 'Select a folder containing files'; 
    answer = questdlg('What do you want to do?', 'Choose source', ...
        sel_files, sel_folder, sel_folder); 
    if strcmpi(answer, sel_files) 
        [file_names, file_path] = uigetfile(fullfile(recent_dir, '*.dat'), 'Select files', 'MultiSelect', 'on');
    else 
        file_path = uigetdir(recent_dir, 'Select folder'); 
        file_names = dir([file_path, '/*.dat']); 
        file_names = {file_names.name}; 
    end
    if ~iscell(file_names), file_names = {file_names}; end 
    file_names = cellfun(@(x) fullfile(file_path, x), file_names, 'uni', 0); 
    recent_dir = file_path; 
else 
    if ~iscell(file_names), file_names = {file_names}; end 
end

cur_tic = tic;
res = cell(length(file_names), 1); 
for file_ind = 1:length(file_names) 
    file_name = file_names{file_ind}; 
    fid = fopen(file_name,'r'); 

    fprintf('Reading %s ', file_name); 
    while ~feof(fid)
         line = fgetl(fid);  
         if ~isempty(regexp(line, '^CONDITION:', 'once'))
             conditions = split(regexprep(line, '^CONDITION:', ''), ';');             
         elseif ~isempty(regexp(line, '^VARIABLES:', 'once'))
             variable_names = split(regexprep(line, '^VARIABLES:', ''));
         elseif ~isempty(regexp(line, '^START_DATA', 'once'))
            matdims = split(fgetl(fid));            
            nrows = str2double(matdims{1});
            ncols = str2double(matdims{2});
            format = [repmat('%f\t', 1, ncols), '%s'];
            data_matrix = textscan(fid, format, nrows);
            data_matrix = double([data_matrix{1:ncols}]);           
         elseif ~isempty(regexp(line, '^END_DATA', 'once'))
             fprintf('... reached end. Current time elapsed = %.2f seconds. \n', toc(cur_tic));
         elseif isempty(line)
             continue; 
         else 
             error('Unexpected line beginnings: %s, at line %d', line, ftell(fid));
         end
    end
    res{file_ind} = struct(...
        'condition', get_condition_struct(conditions), ...
        'data',  get_data_struct(variable_names, data_matrix)); 
    fclose(fid);
end 

res = vertcat(res{:}); 

end

function res = get_condition_struct(conditions)
res = struct; 
for cond_ind = 1:length(conditions) 
    cond = split(conditions{cond_ind}, '=');
    cond_var = cond{1}; 
    cond_val = cond{2}; 
    if ~isnan(str2double(cond_val)), cond_val = str2double(cond_val); end 
    res.(cond_var) = cond_val;
end
end

function res = get_data_struct(variable_names, data_matrix)
if length(variable_names) ~= size(data_matrix, 2) 
    error('Mismatch number of variables and saved data columns');
end

res = struct; 
for var_ind = 1:length(variable_names)
    res.(variable_names{var_ind}) = data_matrix(:, var_ind); 
end
end