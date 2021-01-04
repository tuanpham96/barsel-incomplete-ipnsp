clc; clear; close all; 
run start_up; 

file_prefix = 'full_var-IP-onoff_var-kr-var-abinit_var-norminp_IPnoSP'; 

load(fullfile('data', file_prefix));
fig_path = 'figures/var_2'; 
if ~exist(fig_path, 'dir') 
    mkdir(fig_path);
end

%% Certain prepocessing steps
params_tbl = array2table(param_combs, 'VariableNames', param_labels);

[desired_a0_mat, desired_b0_mat] = meshgrid(unique(params_tbl.a_init), unique(params_tbl.b_init));

%% Define select options for figure
select_opts = struct; 

select_opts.fig = struct('norm_input_opt', 1, 'k_r', 0.1); 

select_opts.col.none_ip = struct('eta_ip_a', 0, 'eta_ip_b', 0);
select_opts.col.gain_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_b', 0);
select_opts.col.bias_ip = struct('eta_ip_a', 0, 'eta_ip_b', max_eta_ip_b);
select_opts.col.full_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_b', max_eta_ip_b);

select_opts.row.complete_train = struct('train_p_inc', 0);
select_opts.row.incomplete_train = struct('train_p_inc', 0.5);

colfields = fieldnames(select_opts.col); 
rowfields = fieldnames(select_opts.row); 


%% Obtain selection 
vis_results = cell(length(rowfields), length(colfields));

for i = 1:length(rowfields)
    for j = 1:length(colfields)
        
        % define selections
        row_fname = rowfields{i};
        col_fname = colfields{j};
        select_conditions = mergefield_struct(select_opts.fig, ...
            select_opts.row.(row_fname), select_opts.col.(col_fname));
        
        % info 
        info_struct = struct;         
        info_struct.select_conditions = select_conditions;
        info_struct.name = [row_fname, '-', col_fname];  
        
        % select
        select_entries = cellfun(@(x) params_tbl.(x) == select_conditions.(x), fieldnames(select_conditions), 'uni', 0); 
        select_entries = all(horzcat(select_entries{:}), 2);
        
        % init 
        init_struct = struct;
        init_struct.a0 = params_tbl.a_init(select_entries);
        init_struct.b0 = params_tbl.b_init(select_entries);
        
        % vis res
        vis_res = get_neccessary_results(train_results(select_entries)); 
        
        % reshape to mat 
        vis_res = structfun(@(v) vec2mat_basedoninits(v, ...
            init_struct.a0, init_struct.b0, desired_a0_mat, desired_b0_mat), ...
            vis_res, 'uni', 0);
        
        init_struct = structfun(@(v) vec2mat_basedoninits(v, ...
            init_struct.a0, init_struct.b0, desired_a0_mat, desired_b0_mat), ...
            init_struct, 'uni', 0);
        
        vis_res = mergefield_struct(vis_res, init_struct);
        
        % save 
        tmp_struct = struct; 
        tmp_struct.info = info_struct;
        tmp_struct.init = init_struct; 
        tmp_struct.res = vis_res;
        
        vis_results{i,j} = tmp_struct; 
        
    end
end

vis_results = cell2mat(vis_results);

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
 
function extracted_results = get_neccessary_results(source_results)
    extracted_results = struct; 
    
    extracted_results.a_f = arrayfun(@(x) x.a.mean(end), source_results);
    extracted_results.a_Df0 = arrayfun(@(x) diff(x.a.mean([1,end])), source_results);

    extracted_results.b_f = arrayfun(@(x) x.b.mean(end), source_results);
    extracted_results.b_Df0 = arrayfun(@(x) diff(x.b.mean([1,end])), source_results);

    extracted_results.dYcomp_f = arrayfun(@(x) x.dY_objvsnoise_test.mean(1,end), source_results);
    extracted_results.dYcomp_Df0 = arrayfun(@(x) diff(x.dY_objvsnoise_test.mean(1,[1,end])), source_results);

    extracted_results.dYincomp_f = arrayfun(@(x) x.dY_objvsnoise_test.mean(2,end), source_results);
    extracted_results.dYincomp_Df0 = arrayfun(@(x) diff(x.dY_objvsnoise_test.mean(2,[1,end])), source_results);
end

function mat = vec2mat_basedoninits(vec, a_v, b_v, a_mat, b_mat)
    mat = zeros(size(a_mat));
    for i = 1:size(mat,1)
        for j = 1:size(mat,2)
            mat(i,j) = vec(a_v == a_mat(i,j) & b_v == b_mat(i,j));
        end
    end
end
