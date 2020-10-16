function [tf, fields_A, fields_B] = similar_fieldlist(A, B, strict) 
fields_A = fieldnames(A); 
fields_B = fieldnames(B); 
tf = isequal(fields_A, fields_B); 

if ~exist('strict', 'var'), strict = true; end

if ~isequal(fields_A, fields_B) && strict
    error('The two structs do not share similar fields');
end

end