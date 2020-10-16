function varargout = return_combination(varargin)
% Examples:
% %% A simple case with `vec` as numeric vectors
%     params_to_vary = {...
%         struct('label', 'N_H', 'vec', [1,50,100]), ...
%         struct('label', 'eta_ip', 'vec', [0, logspace(-5, -1, 5), 0.5, 1]), ...
%         struct('label', 'eta_sp', 'vec', [0, logspace(-5, -1, 5), 0.5, 1]), ...
%         struct('label', 'T', 'vec', [1e3, 1e5]), ...    
%         }; 
%     [combinations, labels, comb_indices] = return_combination(params_to_vary{:}); 
%
% %% A more complicated case, with `vec`'s of mixed types
%     name_of_notA_strvec =  {'OBJ-001', 'OBJ-002', 'NOISE'};
%     genWmethod_strvec = {'all', 'separate'};
%     genWmainfun_structvec = [{struct('type', 'uniform'), struct('type', 'equal')}, ...
%         arrayfun(@(s) struct('type', 'lognorm', 'param', s), [0.5,1,2,5], 'uni', 0)];
% 
%     params_to_vary = {...
%         struct('label', 'name_of_notA', 'vec', {name_of_notA_strvec}), ...
%         struct('label', 'genWmethod', 'vec', {genWmethod_strvec}), ...
%         struct('label', 'genWmainfun', 'vec', {genWmainfun_structvec}), ...
%         struct('label', 'alpha_W', 'vec', [0.6,0.75,0.9,0.95,1]), ...    
%         }; 
% 
%     [combinations, labels, comb_indices] = return_combination(params_to_vary{:}); 


if nargin == 1 && ~isstruct(varargin{1})
    error('Need more than just 1 input that is not a struct');
end

inp_has_param = any(cellfun(@ischar, varargin, 'uni', 1)); 
inp_is_struct = any(cellfun(@isstruct, varargin, 'uni', 1)); 

ip = inputParser; 
default_num = -1; 
valid_num = @(x) isnumeric(x) && length(x) == 1;
ip.addParameter('num', default_num, valid_num);
ip.addParameter('label', {}, @iscell); 

if inp_has_param && inp_is_struct
    error('Cannot be having parameters and a first struct input together');
end

if inp_has_param && ~inp_is_struct 
    ip.parse(varargin{:});     
    vec = {ip.Results.vec}; 
    num = ip.Results.num;
    label = ip.Results.label;

    if num ~= default_num && ~isempty(label) && num ~= length(label)
        err_msg = ['Using one input struct vector and an additional parameter' ...
            ' needs to satisfy either of these following conditions'];
        error('%s\n%s\n%s\n%s\n', err_msg, ...
            '- Has just `num` to tell how many number of repetitions', ...
            '- Has just `label` to tell how the numer of repetitions, output as a struct', ...
            '- Redundantly has both, but their lengths need to be similar');
    end
    
    if num == default_num 
        num = length(label);
    end
end

if inp_is_struct && ~inp_has_param
    vec = cellfun(@(x) x.vec, varargin, 'uni', 0); 
    num = length(vec);
    label = cellfun(@(x) x.label, varargin, 'uni', 0);
end

if ~inp_is_struct && ~inp_has_param
    vec = varargin;
    num = length(vec); 
    label = {}; 
end

    
veclens = cellfun(@(x) 1:length(x), vec, 'uni', 0); 
index_mat = cell(num,1); 
[index_mat{:}] = ndgrid(veclens{:}); 

tmp_cell = cell(num,1);
[tmp_cell{:}] = ndgrid(vec{:}); 
tmp_cell = cellfun(@(x) x(:), tmp_cell, 'uni', 0);

vec_isnumericvec = cellfun(@(v) isnumeric(v) && isvector(v), vec, 'uni', 1);
if all(vec_isnumericvec)
    combo_mat = horzcat(tmp_cell{:});
else
    tmp_cell(vec_isnumericvec) = cellfun(@mat2colcell, tmp_cell(vec_isnumericvec), 'uni', 0);
    combo_mat = horzcat(tmp_cell{:});
end

varargout{1} = combo_mat; 
varargout{2} = label;
varargout{3} = index_mat; 
varargout{4} = struct('combinations', {combo_mat}, 'labels', {label}, 'indices', index_mat);


end