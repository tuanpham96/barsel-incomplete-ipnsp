function s = setfield_structarray(s, field, vals)
vals = num2cell(vals); 
[s.(field)] = vals{:}; 
end