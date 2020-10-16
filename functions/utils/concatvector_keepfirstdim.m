function v = concatvector_keepfirstdim(varargin)
v = varargin{1}; 
if size(v,1) > 1
    all_vectors = cellfun(@to_col_vec, varargin, 'uni', 0); 
    v = vertcat(all_vectors{:}); 
else 
    all_vectors = cellfun(@to_row_vec, varargin, 'uni', 0); 
    v = horzcat(all_vectors{:}); 
end

end 