clc; clear; close all;
run start_up;

graphic_setdefault(25, ...
    'DefaultFigureWindowStyle','normal', ...
    'DefaultAxesTitleFontSize', 1, ...
    'DefaultAxesLabelFontSize', 1.2, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultAxesBox', 'on', ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultFigureWindowStyle','normal');


data_file_prefix = 'full_var-sigmoidwiththres-IP-onoff_var-init_var-psigmoid_IPfasterthanSP';

load(fullfile('data', data_file_prefix));
fig_path = 'figures/var_5';
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

fig_file_prefix = 'sigmoidwiththres-time_at_midT'; 
full_fig_file_prefix = fullfile(fig_path, fig_file_prefix); 

select_time = ceil(size(train_results(1).a.mean,2)/2);
select_time_str = 'T/2';

% select_time = size(train_results(1).a.mean,2);
% select_time_str = 'T';
%% Certain prepocessing steps
params_tbl = array2table(param_combs, 'VariableNames', param_labels);

[desired_a0_mat, desired_theta0_mat] = meshgrid(unique(params_tbl.a_init), unique(params_tbl.theta_init));
desired_init_mat = struct;
desired_init_mat.a0 = desired_a0_mat;
desired_init_mat.theta0 = desired_theta0_mat;

%% Split dY_objvsnoise_test for ease
n_dYtest = size(train_results(1).dY_objvsnoise_test.mean,1);
dYtest_fullstruct = [train_results.dY_objvsnoise_test];
for i = 1:n_dYtest
    struct_arr_ith = arrayfun(@(x) struct('mean', x.mean(i,:), 'sem', x.sem(i,:)), dYtest_fullstruct);
    new_fieldname = sprintf('dYtest%d', i-1);
    for j = 1:length(train_results)
        train_results(j).(new_fieldname) = struct_arr_ith(j); %#ok<SAGROW>
    end
end

%% Define select options for figure
select_opts = struct;

select_opts.col.none_ip = struct('eta_ip_a', 0, 'eta_ip_theta', 0);
select_opts.col.gain_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_theta', 0);
select_opts.col.thres_ip = struct('eta_ip_a', 0, 'eta_ip_theta', max_eta_ip_theta);
select_opts.col.full_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_theta', max_eta_ip_theta);

select_opts.row.train0 = struct('train_p_inc', 0);
select_opts.row.train1 = struct('train_p_inc', 0.2);
select_opts.row.train2 = struct('train_p_inc', 0.4);

colfields = fieldnames(select_opts.col);
rowfields = fieldnames(select_opts.row);

field_list = {'a', 'theta', 'bar_alpha_W', 'dYtest0', 'dYtest1', 'dYtest2'};

%% Latex structs
latex_structs = struct;

latexname_struct = struct;

latexname_struct.none_ip = 'no IP';
latexname_struct.gain_ip = 'gain plast';
latexname_struct.thres_ip = 'threshold plast';
latexname_struct.full_ip = 'full IP';

latexname_struct.train0 = 'p_{\mathrm{inc}}^{\mathrm{train}} = 0';
latexname_struct.train1 = 'p_{\mathrm{inc}}^{\mathrm{train}} = 0.2';
latexname_struct.train2 = 'p_{\mathrm{inc}}^{\mathrm{train}} = 0.4';

latexname_struct.dYtest0 = '\Delta Y | p_{\mathrm{inc}}^{\mathrm{test}} = 0';
latexname_struct.dYtest1 = '\Delta Y | p_{\mathrm{inc}}^{\mathrm{test}} = 0.2';
latexname_struct.dYtest2 = '\Delta Y | p_{\mathrm{inc}}^{\mathrm{test}} = 0.4';

latexname_struct.eta_ip_a = '\eta_a';
latexname_struct.eta_ip_theta = '\eta_{\theta}';
latexname_struct.a = 'a';
latexname_struct.theta = '\theta';
latexname_struct.bar_alpha_W = '\overline{\alpha_W}';

latexfun_struct = struct;
latexfun_struct.text = @(x) x;
latexfun_struct.textbf = @(x) sprintf('\\textbf{%s}',x);
latexfun_struct.textit = @(x) sprintf('\\textit{%s}',x);

latexfun_struct.math = @(x) sprintf('$%s$',x);
latexfun_struct.mathbf = @(x) sprintf('$\\mathbf{%s}$',x);

latexfun_struct.init = @(x) sprintf('$%s_{0}$', x);
latexfun_struct.s = @(x) sprintf('$f = %s$, showing $f(%s)$', x, select_time_str);
latexfun_struct.Ds0 = @(x) sprintf('$f = %s$, showing $f(%s) - f(0)$', x, select_time_str);
latexfun_struct.df = @(x) sprintf('$f = %s$, showing speed $\\langle df_{0 \\rightarrow %s} \\rangle $', x, select_time_str);

latex_structs.names = latexname_struct;
latex_structs.functions = latexfun_struct;

%%
vis_results = plot_selected(full_fig_file_prefix, train_results, params_tbl, field_list, ...
    desired_init_mat, select_opts, select_time, latex_structs);


%% Local functions

function extracted_results = get_neccessary_results(source_results, field_list, select_time)
extracted_results = struct;

for i = 1:length(field_list)
    field_name = field_list{i};
    val_mean_vec = [source_results.(field_name)];
    val_mean_vec = vertcat(val_mean_vec.mean);
    
    extracted_results.([field_name '_s']) = val_mean_vec(:,select_time);
    extracted_results.([field_name '_Ds0']) = val_mean_vec(:,select_time) - val_mean_vec(:,1);
    extracted_results.([field_name '_df']) = mean(diff(val_mean_vec(:,1:select_time),[],2),2);
end
end

function mat = vec2mat_basedoninits(vec, a_v, theta_v, a_mat, theta_mat)
mat = zeros(size(a_mat));
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        mat(i,j) = vec(a_v==a_mat(i,j) & theta_v == theta_mat(i,j));
    end
end
end

function s = rename_with_latex(s, latexname_struct, latexfun_struct)
ind_sep = regexp(s, '_');
if isempty(ind_sep)
    return;
end
ind_sep = ind_sep(end);
var_name = latexname_struct.(s(1:ind_sep-1));
fun_name = latexfun_struct.(s(ind_sep+1:end));
s = fun_name(var_name);
end

function vis_results = plot_selected(full_fig_file_prefix, ...
    train_results, params_tbl, field_list, ...
    desired_init_mat, condition_struct, select_time, latex_structs)
a0_mat = desired_init_mat.a0;
theta0_mat = desired_init_mat.theta0;

latexname_struct = latex_structs.names;
latexfun_struct = latex_structs.functions;

col_conds = condition_struct.col;
row_conds = condition_struct.row;

col_fields = fieldnames(col_conds);
row_fields = fieldnames(row_conds);

ncols = length(col_fields);
nrows = length(row_fields);

vis_results = cell(length(row_fields), length(col_fields));
for i = 1:nrows
    for j = 1:ncols
        row_cond = row_conds.(row_fields{i});
        col_cond = col_conds.(col_fields{j});
        select_conditions = mergefield_struct(row_cond, col_cond);
        
        select_entries = cellfun(@(x) params_tbl.(x) == select_conditions.(x), ...
            fieldnames(select_conditions), 'uni', 0);
        select_entries = all(horzcat(select_entries{:}), 2);
        
        a0_vec = params_tbl.a_init(select_entries);
        theta0_vec = params_tbl.theta_init(select_entries);
        
        vis_res = get_neccessary_results(train_results(select_entries), field_list, select_time);
        
        vis_results{i,j} = structfun(@(v) vec2mat_basedoninits(v, a0_vec, theta0_vec, ...
            a0_mat, theta0_mat), vis_res, 'uni', 0);
        
    end
end

vis_results = cell2mat(vis_results);

bound_structs = arrayfun(@(x) structfun(@(s) get_bound(s), x, 'uni', 0), vis_results);
bound_structs = structfun(@(x) vertcat(x{:}), structarray_to_struct(bound_structs, 0), 'uni', 0);
bound_structs = structfun(@(x) [min(x(:,1)), max(x(:,2))], bound_structs, 'uni', 0);

fullfield_list = fieldnames(vis_results);
cmap = return_colorbrewer('RdBu_r', 100);

rename_fun_short = @(s) rename_with_latex(s, latexname_struct, latexfun_struct);

% x_vec = unique(params_tbl.a_init);
% y_vec = unique(params_tbl.theta_init);

x_lbl = rename_fun_short('a_init');
y_lbl = rename_fun_short('theta_init');

t0 = tic; 
progress_strnchar = fprintf('Starting plotting ...'); 

nfigs = length(fullfield_list);
for k = 1:nfigs
    field_name = fullfield_list{k};
    fig_filename = sprintf('%s-%s', full_fig_file_prefix, field_name);
    
    shared_clim_bound = bound_structs.(field_name);
    figure('units', 'normalized', 'position', [0.1,0.1,0.7,0.9]);
    ax_pos = tight_subplot(nrows,ncols,[.03 .03],[0.07 0.1],[0.08 0.1], false);
    
    fig_description = rename_fun_short(field_name);
    annotation('textbox', 'units', 'normalized', 'position', [0,0.92,1,0.08], ...
        'string', fig_description, 'linestyle', 'none', ...
        'fontsize', 30, 'interpreter', 'latex', ...
        'horizontalalignment', 'center', 'verticalalignment', 'top');
    for i = 1:nrows
        for j = 1:ncols
            
            row_field = row_fields{i};
            col_field = col_fields{j};
            
            axes('units', 'normalized', 'position', ax_pos{i,j}, ...
                'xtick', '', 'ytick', '');
            hold 'on';
            image_with_strict_limits(vis_results(i,j).(field_name)');
            
            if i == 1, title(rename_fun_short([col_field, '_textbf'])); end
            
            if j == 1, ylabel({rename_fun_short([row_field, '_mathbf']),y_lbl}); end
            
            if i == nrows, xlabel(x_lbl), end
            
            if i == nrows && j == ncols
                % set(gca, 'xtick', [1, length(x_vec)], 'xticklabels', x_vec([1,end]), ...
                %    'ytick', [1, length(y_vec)], 'yticklabels', y_vec([1,end]));
                cbar = colorbar;
                cbar.Position = cbar.Position .* [1,1,1,0.8] + [0.06,0.01,0,0];
                cbar.FontSize = 20;
            end
        end
    end
    
    colormap(cmap);
    ax_list = findall(gcf, 'type', 'axes');
    arrayfun(@(ax) daspect(ax, [1,1,1]), ax_list);
    share_caxis(ax_list, shared_clim_bound);
    
    export_fig(fig_filename, '-r300', '-p0.02');
    close;
    
    t1 = toc(t0)/60; 
    fprintf(repmat('\b', 1, progress_strnchar));
    progress_strnchar = fprintf('(%d/%d) - field "%s" plotted and saved. Elapsed %.1f min. ETA %.1f min.', ...
        k, nfigs, field_name, t1, (nfigs/k - 1)*t1);
    
end

fprintf('\nDone\n');
end
