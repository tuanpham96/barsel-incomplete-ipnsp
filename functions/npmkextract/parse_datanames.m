function file_list = parse_datanames(folder, file_types)
if ~exist('file_types', 'var')
    file_types = {'*.mat', '*.bin', '*.npy', '*.dat', '*.mda'};  % add more along the way if needed
end 

file_list = cellfun(@(x) dir(fullfile(folder, x)), file_types, 'uni', 0); 
file_list = arrayfun(@(x) x.name, vertcat(file_list{:}), 'uni', 0); 
file_list = cellfun(@(x) parse_onename(x, folder), file_list, 'uni', 1); 
end

function db = parse_onename(file_name, folder) 
file_name_suffix_split = split(file_name, '.'); 
file_name_split = split(file_name_suffix_split{1}, '_'); 

db = struct; 
db.type = file_name_suffix_split{2}; 

if length(file_name_split) < 4 || length(file_name_split) > 5
    error(['The file name "%s" needs to follow this format: \n' ...
        '<subject>_<condition>_<condition_index>_<analysis>(_<analysis_index>).<type>'], file_name);
end 
db.subject = file_name_split{1}; 
db.condition = [file_name_split{2} file_name_split{3}];

if length(file_name_split) == 4
    analysis = file_name_split{4};
else 
    analysis = [file_name_split{4} file_name_split{5}]; 
end 

db.analysis = analysis; 
db.path = fullfile(folder, file_name); 
db.identifier = file_name_suffix_split{1}; 
end