function res = ifelse(condition, true_value, varargin)
% If-else shorthand, borrowed ideas from many other languages 
% also, there is a function available in MATLAB Exchange: 
%   https://www.mathworks.com/matlabcentral/fileexchange/28587-if-then-else-shortcut
% However, this one does not consider the returned values to be expressions
% to be evaluated by 'eval*'-like functions to be debug-safe-function
% used as: 
%   ifelse(condition, true_value, false_value)
%   ifelse(condition_1, true_value_1, condition_2, true_value_2, false_value) 
%   ifelse([<condition_i, true_value_i>,...], false_value) to string along
%       multiple ifelse

if mod(length(varargin),2) ~= 1
    error('Needs pairs of (condition, true_value) and followed by the last false_value');
end
if length(condition) ~= 1, error('Only scalar boolean is allowed'); end 

if condition 
    res = true_value; 
else
    if length(varargin) == 1
        res = varargin{1};
    else 
        res = ifelse(varargin{:});
    end
end