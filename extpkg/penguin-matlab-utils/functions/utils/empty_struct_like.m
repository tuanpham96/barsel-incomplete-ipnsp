function structobj = empty_struct_like(template_with_fields, init_value)
if isstruct(template_with_fields)
    fields = fieldnames(template_with_fields)';
elseif iscell(template_with_fields) 
    fields = unique(to_row_vec(template_with_fields));
else 
    error('Only accept either a struct or a cell list of fields');
end

empties = cell(1,length(fields)); 
if exist('init_value', 'var'), empties = cellfun(@(x) init_value, empties, 'uni', 0); end
cat_field_empty = [fields; empties]; 
structobj = struct(cat_field_empty{:}); 
end