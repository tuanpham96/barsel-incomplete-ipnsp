function c = mat2colcell(m)
% Matrix to colum cell 
c = mat2cell(m, ones(size(m,1),1)); 
end