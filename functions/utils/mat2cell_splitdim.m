function c = mat2cell_splitdim(m,d)
% Matrix to row cell 
sz = num2cell(size(m));
sz{d} = ones(sz{d}, 1); 
c = squeeze(cellfun(@squeeze, mat2cell(m, sz{:}), 'uni', 0)); 
end