function save_figure(varargin)
first_inp = varargin{1}; 
curr_argind = 2; 
if strcmp(first_inp.Type, 'figure')
    fig = first_inp; 
else
    fig = gcf; 
    curr_argind = 1; 
end 
figr_name = varargin{curr_argind}; curr_argind = curr_argind + 1;  
figr_path = varargin{curr_argind}; curr_argind = curr_argind + 1;  

save_figure_with_param(fig, figr_name, figr_path, varargin{curr_argind:end}); 

end

function save_figure_with_param(fig, figr_name, figr_path, varargin)

fprintf('Saving %s. Please wait ... ', figr_name);

if isempty(varargin) 
    print_opts = {'-dpdf', '-painters'}; 
else
    print_opts = varargin;
end

print(fig, fullfile(figr_path, figr_name), print_opts{:});
fprintf('Done saved.\n');

if any(contains(print_opts, 'pdf', 'IgnoreCase', true))
    fprintf('Opening ...\n');
    open([fullfile(figr_path, figr_name) '.pdf']);
end

end