clc; clear; close all; 
run start_up; 

file_prefix = 'full_var-sigmoidwiththres-IP-onoff_var-init_var-psigmoid_IPfasterthanSP'; 

load(fullfile('data', file_prefix));
fig_path = 'figures/var_5'; 
if ~exist(fig_path, 'dir') 
    mkdir(fig_path);
end

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
% select_time = ceil(size(train_results(1).a.mean,2)/2);
% select_time_str = 'T/2'; 

select_time = size(train_results(1).a.mean,2);
select_time_str = 'T'; 

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

field_list = {'a', 'theta', 'bar_alpha_W', 'dYtest0', 'dYtest2'}; 

latexname_struct = struct; 

latexname_struct.none_ip = '\mathrm{no IP}';
latexname_struct.gain_ip = '\mathrm{gain plast}';
latexname_struct.thres_ip = '\mathrm{threshold plast}';
latexname_struct.full_ip = '\mathrm{full IP}';

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
latexname_struct.bar_alpha_W = '\overline_{\alpha_W}';


latexfun_struct = struct; 
latexfun_struct.s = @(x) sprintf('$f = %s$, showing $f(%s)$', x, select_time_str);
latexfun_struct.Ds0 = @(x) sprintf('$f = %s$, showing $f(%s) - f(0)$', x, select_time_str);
latexfun_struct.df = @(x) sprintf('$f = %s$, showing speed $\\langle df_{0\\rightarrow%s}\\rangle$', x, select_time_str);

return;
vis_results = plot_selected(train_results, params_tbl, field_list, desired_init_mat, select_opts, select_time);
return;
%% Plot
nrows = length(rowfields);
ncols = length(colfields);

plt_field_prefix = {'a', 'b', 'dYcomp', 'dYincomp'};
plt_field_suffices = {'f', 'Df0'};

for i_pltf_suf = 1:length(plt_field_suffices)
plt_field_suffix = plt_field_suffices{i_pltf_suf};

cnt_splt = 1; 
outer_boxcolor = 'none';

if strcmpi(plt_field_suffix, 'f') 
    pseudo_str_pltsuff = 'final values';
    fig_name_pltsuff = 'f'; 
else    
    pseudo_str_pltsuff = 'final - init';
    fig_name_pltsuff = 'Df0'; 
end

if select_opts.fig.norm_input_opt
    pseudo_str_norminp = 'with norm';
    fig_name_norminp = 'norminp';
else    
    pseudo_str_norminp = 'without norm';
    fig_name_norminp = 'nonorminp';
end

if select_opts.fig.k_r == 0.1
    pseudo_str_selectivity = 'high';
    fig_name_selectivity = 'highsel';
else    
    pseudo_str_selectivity = 'medium';
    fig_name_selectivity = 'medsel';
end

fig_description = sprintf('plot %s of variables, %s input, %s selectivity', ...
    pseudo_str_pltsuff, pseudo_str_norminp, pseudo_str_selectivity);
fig_name = sprintf('%s-plot-%s-var-%s_%s', file_prefix, ...
    fig_name_pltsuff, fig_name_norminp, fig_name_selectivity);

figure;

annotation('textbox', 'units', 'normalized', 'position', [0,0.92,1,0.08], ...
    'string', fig_description, 'linestyle', 'none', ...
    'fontsize', 25, 'interpreter', 'latex', ...
    'horizontalalignment', 'center', 'verticalalignment', 'top'); 
for i = 1:length(rowfields)
    for j = 1:length(colfields)
        
        row_fname = rowfields{i};
        col_fname = colfields{j};
        
        main_ax = subplot(nrows, ncols, cnt_splt); 
        hold on; 
        
        title(regexprep([row_fname, '-', col_fname], '_', '-'));
        main_ax_pos = get(main_ax, 'position'); 
        
        subax_commonpos = main_ax_pos .* [1,1,0.4,0.4]; 
        [subax_horzshift, subax_vertshift]= meshgrid([0,0.10], [0.15,-0.03]); 
        
        for k = 1:length(plt_field_prefix)
            sub_ax_pos = subax_commonpos + [subax_horzshift(k),subax_vertshift(k),0,0]; 
            sub_ax_handle = axes(gcf, 'Units', 'normalized', 'Position', sub_ax_pos, 'Tag', plt_field_prefix{k});
            hold on;
            
            mat_toplot = vis_results(i,j).res.([plt_field_prefix{k} '_' plt_field_suffix]); 
            image_with_strict_limits(sub_ax_handle, mat_toplot');
            set(sub_ax_handle, 'xtick', '', 'ytick', '', 'box', 'on', 'linewidth', 3); 
            
            title(plt_field_prefix{k});
            
            if cnt_splt == 1
                
                if k == 1
                    xlabel('$a_0$');
                    ylabel('$b_0$');
                end
                
                cbar = colorbar;
                cbar.Position = cbar.Position .* [1,1,0.8,0.5] + [0.05,0.01,0,0]; 
                cbar.FontSize = 12;
            end
        end
        
        set(main_ax, 'Position', main_ax_pos .* [1,1,1.2,1.2] + [-0.015,-0.07,0,0], ...
            'color', 'none', 'xtick', '', 'ytick', '', 'box', 'on', 'linewidth', 1.5, ...
            'xcolor', outer_boxcolor, 'ycolor', outer_boxcolor);
        cnt_splt = cnt_splt + 1; 
        
    end
end

arrayfun(@(ax) daspect(ax, [1,1,1]), findall(gcf, 'type', 'axes'));

all_colorlims = cellfun(@(tagname) ...
    cell2mat(arrayfun(@(ax) caxis(ax), findall(gcf, 'type', 'axes', 'tag', tagname), 'uni', 0)), ...
    plt_field_prefix, 'uni', 0);
all_colorlims = cellfun(@(x) [min(x(:,1)), max(x(:,2))], all_colorlims, 'uni', 0);

if strcmpi(plt_field_suffix, 'Df0')
    all_colorlims = cellfun(@(x) [-1,1] * max(abs(x)), all_colorlims, 'uni', 0);
    colormap(return_colorbrewer('RdBu_r', 500));
else
    all_colorlims = cellfun(@(x) x+[-1,1]*0.01*diff(x), all_colorlims, 'uni', 0);
    colormap(return_colorbrewer('YlGnBu', 100)*0.98);
end

cellfun(@(tagname, clim) ...
    arrayfun(@(ax) caxis(ax, clim), findall(gcf, 'type', 'axes', 'tag', tagname)), ...
    plt_field_prefix, all_colorlims);
        
% export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
% close; 
end

%% Local functions
 
function extracted_results = get_neccessary_results(source_results, field_list, select_time)
    extracted_results = struct; 
    
    for i = 1:length(field_list)
        field_name = field_list{i};
        val_mean_vec = [source_results.(field_name)]; 
        val_mean_vec = vertcat(val_mean_vec.mean);
        
%         extracted_results.([field_name '_s']) = val_mean_vec(:,select_time); 
        extracted_results.([field_name '_Ds0']) = val_mean_vec(:,select_time) - val_mean_vec(:,1); 
%         extracted_results.([field_name '_df']) = mean(diff(val_mean_vec(:,1:select_time),[],2),2); 
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

function vis_results = plot_selected(train_results, params_tbl, field_list, desired_init_mat, condition_struct, select_time)
a0_mat = desired_init_mat.a0;
theta0_mat = desired_init_mat.theta0; 

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

fullfield_list = fieldnames(vis_results);
cmap = return_colorbrewer('RdBu_r', 100); 

for k = 1:length(fullfield_list)
    field_name = fullfield_list{k};
    
    figure; cnt = 1;
    for i = 1:nrows
        for j = 1:ncols
            subplot(nrows, ncols, cnt); hold on;
            image_with_strict_limits(vis_results(i,j).(field_name)');
            cnt = cnt + 1;
            title(regexprep([field_name ' - ' col_fields{j}], '_', '-'));
            ylabel(regexprep(row_fields{i}, '_', '-'));
            
            set(gca, 'xtick', '', 'ytick', '');
            colorbar; 
        end
    end
    
    colormap(cmap);
    arrayfun(@(ax) daspect(ax, [1,1,1]), findall(gcf, 'type', 'axes'));

end
end
