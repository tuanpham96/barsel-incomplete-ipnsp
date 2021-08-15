clc; clear; close all;
run start_up;

graphic_setdefault(25, ...
    'DefaultFigureWindowStyle','normal', ...
    'DefaultAxesTitleFontSize', 1.1, ...
    'DefaultAxesLabelFontSize', 1.1, ...
    'DefaultAxesLineWidth', 2, ...
    'DefaultTextInterpreter', 'latex', ...
    'DefaultLegendInterpreter', 'latex', ...
    'DefaultFigureWindowStyle','normal');


data_file_prefix = 'full_var-sigmoidwiththres-IP-onoff_var-init_var-psigmoid_IPfasterthanSP';

load(fullfile('data', data_file_prefix));
fig_path = 'figures/var_5';
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%% Certain prepocessing steps
params_tbl = array2table(param_combs, 'VariableNames', param_labels);
unq_a0 = unique(params_tbl.a_init);
unq_theta0 = unique(params_tbl.theta_init);
[desired_a0_mat, desired_theta0_mat] = meshgrid(unq_a0, unq_theta0);
desired_init_mat = struct;
desired_init_mat.a0 = desired_a0_mat;
desired_init_mat.theta0 = desired_theta0_mat;

%%
% eql_or_not = false(length(train_results),1);
% for i = 1:length(train_results)
%     prm = table2struct(params_tbl(i,:)); 
%     res = train_results(i).opts;
%     eql_or_not(i) = all(cellfun(@(x) res.(x) == prm.(x), fieldnames(prm), 'uni', 1));
% end

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

%% Selecting (a0,theta0) pairs
ind_theta0_vec = [3, 8, 7];
ind_a0_vec = [2, 3, 8]; 

sel_color_inds = [3,5,8];
cmaps = cellfun(@(x) return_colorbrewer(x), {'Greys', 'PuRd', 'GnBu'}, 'uni',0);
cmaps = cellfun(@(x) x(sel_color_inds,:), cmaps, 'uni', 0);

linstyles = {':', '-', '--'};

selected_fields = {'a','theta','bar_alpha_W','dYtest0','dYtest1','dYtest2'}; 
removed_fields = setdiff(fieldnames(train_results), selected_fields);

ip_cond_struct = struct; 
ip_cond_struct.none_ip = struct('eta_ip_a', 0, 'eta_ip_theta', 0);
ip_cond_struct.gain_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_theta', 0);
ip_cond_struct.thres_ip = struct('eta_ip_a', 0, 'eta_ip_theta', max_eta_ip_theta);
ip_cond_struct.full_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_theta', max_eta_ip_theta);
ip_conds = fieldnames(ip_cond_struct); 

train_cond_struct = struct;
train_cond_struct.train0 = struct('train_p_inc', 0);
train_cond_struct.train1 = struct('train_p_inc', 0.2);
train_cond_struct.train2 = struct('train_p_inc', 0.4);
train_conds = fieldnames(train_cond_struct);

num_init_select = length(ind_a0_vec);
selected_results = cell(num_init_select,1);

for i = 1:num_init_select
    a0 = unq_a0(ind_a0_vec(i)); 
    theta0 = unq_theta0(ind_theta0_vec(i));
    
    common_select_conditions = struct; 
    common_select_conditions.a_init = a0; 
    common_select_conditions.theta_init = theta0;
    
    tmp_res = struct(); 
    
    for j1 = 1:length(ip_conds)
        ip_cond = ip_conds{j1};
        
        tmp_res_j1 = struct; 
        for j2 = 1:length(train_conds)
            train_cond = train_conds{j2};
            
            select_conditions = mergefield_struct(common_select_conditions, ...
                ip_cond_struct.(ip_cond), train_cond_struct.(train_cond));
            
            select_entry = cellfun(@(x) params_tbl.(x) == select_conditions.(x), ...
                fieldnames(select_conditions), 'uni', 0);            
            select_entry = all(horzcat(select_entry{:}), 2);
            tmp_res_j1.(train_cond) = rmfield(train_results(select_entry),removed_fields);
        end
        
        tmp_res_j2 = struct; 
        for j2 = 1:length(selected_fields)
            sel_field = selected_fields{j2};
            tmp_j2 = structfun(@(x) x.(sel_field), tmp_res_j1, 'uni', 0);
            tmp_j2 = cellfun(@(x) tmp_j2.(x).mean, train_conds, 'uni', 0);
            tmp_res_j2.(sel_field) = vertcat(tmp_j2{:});
        end
        
        tmp_res.(ip_cond) = tmp_res_j2;
    end
    
    selected_results{i} = tmp_res;
end

%% Latex structs
latex_structs = struct;

latexname_struct = struct;

latexname_struct.none_ip = 'no IP';
latexname_struct.gain_ip = 'gain plast';
latexname_struct.thres_ip = 'threshold plast';
latexname_struct.full_ip = 'full IP';

latexname_struct.train = 'p_{\mathrm{inc}}^{\mathrm{train}}';
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
latexfun_struct.s = @(x) sprintf('$%s$ at $%s$', x, select_time_str);
latexfun_struct.Ds0 = @(x) sprintf('$f = [%s]$, showing $f(%s) - f(0)$', x, select_time_str);
latexfun_struct.df = @(x) sprintf('$f = [%s]$, showing speed $\\langle df_{0 \\rightarrow %s} \\rangle $', x, select_time_str);

latex_structs.names = latexname_struct;
latex_structs.functions = latexfun_struct;

%%
nrows = 3; 
ncols = length(ip_conds); 
select_field_figs = {selected_fields(1:3), selected_fields(4:end)}; 
fignames = {'selectvis-net', 'selectvis-dY'};

for i = 1:length(select_field_figs)
    sel_field_list = select_field_figs{i}; 
    
    figure;
    ax_pos = tight_subplot(nrows,ncols,[.03 .03],[0.1 0.08],[0.08 0.1], false);
    [ax_cols,ax_rows] = meshgrid(1:ncols,1:nrows);
    [fields_1, fields_2] =  meshgrid(ip_conds,sel_field_list); 
    for j = 1:numel(fields_1)

        f1 = fields_1{j};
        f2 = fields_2{j};
        
        ax_cmap = vertcat(cmaps{:});
        axes('units', 'normalized', 'position', ax_pos{j}, ...
            'colororder', ax_cmap, ...
            'Tag', sprintf('%s-%s', f1, f2), 'box', 'off');
        hold on; 
        plt_mat = cellfun(@(x) x.(f1).(f2)', selected_results, 'uni', 0);
        plt_mat = horzcat(plt_mat{:});        
        plot(plt_mat, 'linewidth', 3);
        
        xline((def_opts.num_train * def_opts.p_train_complete)/def_opts.subsampled, ':k', 'linewidth', 2);
        
        row_ind = ax_rows(j);
        col_ind = ax_cols(j);
        if row_ind == 1
            title(latexfun_struct.textbf(latex_structs.names.(f1)));
        end
        if col_ind == 1
            ylabel(sprintf('$%s$', latex_structs.names.(f2)));
        else
            set(gca, 'ycolor', 'none');
        end
        if row_ind == nrows
            xlabel(sprintf('\\# steps $\\times %d$', def_opts.subsampled)); 
        else
            set(gca, 'xcolor', 'none');
        end
        
    end
    
    
    linkaxes(findall(gcf, 'type', 'axes'), 'x'); 
    for j = 1:length(sel_field_list)
        field_name = sel_field_list{j};
        ax_linky = findall(gcf, 'type', 'axes');
        ind_ax = regexp(get(ax_linky,'tag'), sprintf('-%s$', field_name), 'tokens');
        ax_linky = ax_linky(cellfun(@(x) ~isempty(x), ind_ax, 'uni', 1));
        linkaxes(ax_linky, 'y');
    end
    despline('all');
    
    fig_name = fullfile(fig_path, fignames{i});
    
    export_fig(fig_name, '-r300', '-p0.02');
    close;
end

%%
ncols = 8;

figure(...
    'Units', 'normalized', ...
    'Color', 'none',...
    'Position', [0.1,0.1,0.6,0.5], ...
    'DefaultAxesFontSize', 35, ...
    'DefaultAxesTitleFontSize', 1.1, ...
    'DefaultAxesLabelFontSize', 1.6, ...
    'DefaultAxesLineWidth', 2);
ax = subplot(1,ncols+1,1:(ncols-num_init_select)); hold on; 
init_mat = zeros(length(unq_a0), length(unq_theta0));
cmap_init = cellfun(@(x) x(end,:), cmaps, 'uni', 0);
cmap_init = [0.95,0.95,0.95;vertcat(cmap_init{:})];
for i = 1:num_init_select
    init_mat(ind_a0_vec(i), ind_theta0_vec(i)) = i;
end
image(unq_a0, unq_theta0, init_mat'); hold on; 
colormap(gca, cmap_init); pbaspect([1,1,1]);
xlim(unq_a0([1,end]) + 0.5*[-1,1]'*diff(unq_a0([1,2]))); 
ylim(unq_theta0([1,end]) + 0.5*[-1,1]'*diff(unq_theta0([1,2]))); 
xlabel(latexfun_struct.init(latex_structs.names.a));
ylabel(latexfun_struct.init(latex_structs.names.theta));
ax.Position = ax.Position.*[1,1,0.7,0.7] + [0.05,0.2,0,0];
despline;

num_cbar_ticks = length(train_conds); 
cbar_ticklbls = unique(params_tbl.train_p_inc);
for i = 1:num_init_select
    ax = subplot(1,ncols+1,ncols-num_init_select+i); 
    set(ax, 'visible', 'off');
    cbar = colorbar(ax); 
    colormap(ax, cmaps{i});
    cbar.Position = cbar.Position .* [1,1,8,0.7];
    cbar.FontSize = 40;
    cbar.Box = 'on';
    cbar.LineWidth = 3;
    if i == num_init_select
        cbar_ticks = linspace(0,1,num_cbar_ticks+1);
        cbar_ticks = 0.5*(cbar_ticks(1:end-1) + cbar_ticks(2:end)); 
        cbar.Ticks = cbar_ticks;
        cbar.TickLabels = cbar_ticklbls;        
    else
        cbar.Ticks = '';
    end
    if i == 2
        title(cbar, sprintf('$%s$', latexname_struct.train), 'interpreter', 'latex', 'fontsize', 50);
    end
    
end
fig_name = fullfile(fig_path, 'selectvis-sup');
export_fig(fig_name, '-r300', '-p0.02');


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
