function C = minus(A, B)
[~, field_list, ~] = similar_fieldlist(A, B, true); 
C = struct(A); 
for i = 1:length(field_list)
    field_name = field_list{i}; 
    res = getfield_structarray(A, field_name) - getfield_structarray(B, field_name); 
    C = setfield_structarray(C, field_name, res); 
end 

end