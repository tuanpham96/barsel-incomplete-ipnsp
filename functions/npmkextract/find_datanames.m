function filtered_filelist = find_datanames(file_list, varargin)
if ischar(file_list) 
    file_list = parse_datanames(file_list); 
end 

nvar = length(varargin); 
if mod(nvar, 2) ~= 0 || nvar < 2 
    error('Need pairs of additional inputs to find data');
end 
if ~all(isfield(file_list, varargin(1:2:end)))
    error('Not all the variable names are available in `file_list` struct');
end 

nvar = nvar/2; % actual variables to find; 

filtered_filelist = file_list; 
for i = 1:nvar 
    var_name = varargin{2*i-1};
    var_pattern = varargin{2*i}; 
    
%     ind_contains = contains({filtered_filelist.(var_name)}, var_pattern);
    ind_contains = cellfun(@(x) ~isempty(x), ...
        regexp({filtered_filelist.(var_name)}, var_pattern), 'uni', 1); 
    filtered_filelist = filtered_filelist(ind_contains);
end 

end