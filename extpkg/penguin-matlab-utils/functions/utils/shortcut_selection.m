function res = shortcut_selection(list, ttl, return_type, select_mode)
% shortcut_selection: to prompt a selection from a list using `listdialog`
%   list:           a cell list of options
%   ttl:            title of prompt 
%   return_type:    [optional, default='VAL'] 
%                   'VAL':  return value
%                   'IND':  return index 
%   select_mode:    [optional, default='SINGLE']
%                   'SINGLE':   single selection 
%                   'MUL':      multiple selections
if ~exist('return_type', 'var'), return_type = 'VAL'; end
return_type = regexprep(return_type, {'return','_',' '}, ''); 
switch upper(return_type) 
    case {'IND', 'INDEX', 'INDICES', 'IDX'}, return_type = 'IND'; 
    case {'VAL', 'VALS', 'VALUE', 'VALUES'}, return_type = 'VAL'; 
    otherwise, error('Only return_type=''IND'' or ''VAL'' are accepted, not %s', return_type); 
end

if ~exist('select_mode', 'var'), select_mode = 'SINGLE'; end
select_mode = regexprep(select_mode, {'selection','select','_',' '}, ''); 
switch upper(select_mode) 
    case 'SINGLE', select_mode = 'SINGLE'; 
    case {'MULTIPLE', 'MUL'}, select_mode = 'MULTIPLE'; 
    otherwise, error('Only return_type=''single/multiple'' are accepted, not %s', select_mode); 
end

ind = listdlg('ListString', list, 'ListSize', [300, 300], ...
    'SelectionMode', lower(select_mode), 'PromptString', ttl);

if strcmpi(return_type, 'IND')
    res = ind; 
else
    if isempty(ind), second_cond_value = []; else, second_cond_value = list{ind}; end
    res = ifelse(isempty(ind), [], ... % return empty if nothing is selected 
        length(ind) == 1 && strcmpi(select_mode, 'single'), second_cond_value, ... % return the singly-selected string 
        list(ind)); % if select_multiple + something is selected 
end

end