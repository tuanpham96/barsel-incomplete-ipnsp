function alpha = adjust_alpha(alpha, num_nulls, multiadj_meth)
% multicomparison adjustment of alpha 
% TODO: add more 
switch upper(multiadj_meth) 
    case 'NONE' 
    case {'BONFERRONI', 'BONF' } 
        alpha = alpha/num_nulls;
    otherwise 
        error('The multi-comparison method %s is not accepted', multiadj_meth);
end

end