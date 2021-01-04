clc; clear; close all; 
run start_up; 

file_prefix = 'full_var-IP-onoff_var-kr-var-abinit_var-norminp_IPnoSP'; 

global params_tbl train_results select_opts colfields linefields rowfields max_ncolors
global ctrl_fig plt_fig desired_a0_mat desired_b0_mat glob_cmap
global rng_of_a rng_of_b
load(fullfile('data', file_prefix));

max_ncolors = 9;
glob_cmap = return_colorbrewer('Set1',max_ncolors); 
%% Certain prepocessing steps
params_tbl = array2table(param_combs, 'VariableNames', param_labels);


[desired_a0_mat, desired_b0_mat] = meshgrid(unique(params_tbl.a_init), unique(params_tbl.b_init));

rng_of_a = arrayfun(@(x) [min(x.a.mean),max(x.a.mean)], train_results, 'uni', 0);
rng_of_a = vertcat(rng_of_a{:});
rng_of_a = [min(rng_of_a(:,1)), max(rng_of_a(:,2))];

rng_of_b = arrayfun(@(x) [min(x.b.mean),max(x.b.mean)], train_results, 'uni', 0);
rng_of_b = vertcat(rng_of_b{:});
rng_of_b = [min(rng_of_b(:,1)), max(rng_of_b(:,2))];

%% Define select options for figure
select_opts = struct; 

select_opts.fig = struct('norm_input_opt', 1, 'k_r', 0.4, ...
                            'a_init', nan, ...
                            'b_init', nan); 

select_opts.col.none_ip = struct('eta_ip_a', 0, 'eta_ip_b', 0);
select_opts.col.gain_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_b', 0);
select_opts.col.bias_ip = struct('eta_ip_a', 0, 'eta_ip_b', max_eta_ip_b);
select_opts.col.full_ip = struct('eta_ip_a', max_eta_ip_a, 'eta_ip_b', max_eta_ip_b);

select_opts.line.complete_train = struct('train_p_inc', 0);
select_opts.line.incomplete_train = struct('train_p_inc', 0.5);

colfields = fieldnames(select_opts.col); 
linefields = fieldnames(select_opts.line); 
rowfields = {'a', 'b', 'dYcomp', 'dYincomp'};

%% Plot 
ctrl_fig = figure('units', 'normalized', 'position', [0.8,0.1,0.2,0.3]); 
plt_fig = figure('units', 'normalized', 'position', [0,0,0.9,1]);

figure(ctrl_fig);
colormap(glob_cmap); 
mat = nan(size(desired_a0_mat)); 
image(mat', 'AlphaData', ~isnan(mat)');
style_ctrl_fig;
set(gcf, 'KeyPressFcn', @(~,~) choose_pixel);

%% Local functions

function vis_results = obtain_selection(a_init, b_init)

    global params_tbl train_results select_opts colfields linefields;
    select_opts.fig.a_init = a_init; 
    select_opts.fig.b_init = b_init;
    
    vis_results = cell(length(colfields), length(linefields));

    for i = 1:length(linefields)
        for j = 1:length(colfields)

            % define selections
            line_fname = linefields{i};
            col_fname = colfields{j};
            select_conditions = mergefield_struct(select_opts.fig, ...
                select_opts.line.(line_fname), select_opts.col.(col_fname));

            % info
            info_struct = struct;
            info_struct.select_conditions = select_conditions;
            info_struct.name = [line_fname, '-', col_fname];

            % select
            select_entries = cellfun(@(x) params_tbl.(x) == select_conditions.(x), fieldnames(select_conditions), 'uni', 0);
            select_entries = all(horzcat(select_entries{:}), 2);

            % init
            init_struct = struct;
            init_struct.a0 = params_tbl.a_init(select_entries);
            init_struct.b0 = params_tbl.b_init(select_entries);

            % vis res
            vis_res = get_neccessary_results(train_results(select_entries));

            % save
            tmp_struct = struct;
            tmp_struct.info = info_struct;
            tmp_struct.init = init_struct;
            tmp_struct.res = vis_res;

            vis_results{j,i} = tmp_struct;

        end
    end

    vis_results = cell2mat(vis_results);
end

function extracted_results = get_neccessary_results(source_results)
    extracted_results = struct; 
    
    extracted_results.a = source_results.a.mean;
    extracted_results.b = source_results.b.mean;

    extracted_results.dYcomp = source_results.dY_objvsnoise_test.mean(1,:);
    extracted_results.dYincomp = source_results.dY_objvsnoise_test.mean(2,:);
end


function choose_pixel
    global ctrl_fig desired_a0_mat desired_b0_mat glob_cmap
    if ~strcmpi(get(ctrl_fig, 'currentcharacter'), 'c')
        return; 
    end
    
    figure(ctrl_fig);
    
    ax = findall(ctrl_fig, 'type', 'axes');
    mat = get(ax, 'children').CData';
    
    [x_coord,y_coord] = ginput(1); 
    x_coord = min(max(round(x_coord), 1), size(mat,1));
    y_coord = min(max(round(y_coord), 1), size(mat,2));
    
    mat_val = mat(x_coord, y_coord);
    
    max_mat = 0; 
    if any(~isnan(mat(:)))
        max_mat = max(mat(:), [], 'omitnan');
    end
    
    if isnan(mat_val) 
        mat_val = max_mat+1;
        line_color = glob_cmap(mat_val,:);
    else
        mat_val = nan;  
        line_color = nan; 
    end
    
    a_init = desired_a0_mat(y_coord, x_coord);
    b_init = desired_b0_mat(y_coord, x_coord);
    
    mat(x_coord, y_coord) = mat_val;
    
    cla; hold on; 
    image(ax, mat', 'AlphaData', ~isnan(mat'));
    style_ctrl_fig;

    plot_fig_after_choose(a_init, b_init, line_color); 
    
    figure(ctrl_fig);
end

function style_ctrl_fig
    global ctrl_fig glob_cmap 
    figure(ctrl_fig);
    caxis([1,size(glob_cmap,1)]);
    set(gca, 'color', [1,1,1]*0.95, 'ydir', 'normal')
    box on; daspect([1,1,1]);
    xlabel('$a_0$');
    ylabel('$b_0$');
    title('ctrl fig');
    
end

function plot_fig_after_choose(a_init, b_init, line_color)
    global colfields rowfields 
    global plt_fig 
    global rng_of_a rng_of_b

    vis_results = obtain_selection(a_init, b_init);

    figure(plt_fig);
    cnt_splt = 1; 

    nrows = length(rowfields);
    ncols = length(colfields);

    
    tag_init = sprintf('a=%.3f;b=%.4f', a_init, b_init);
    if all(isnan(line_color))
        delete(findall(gcf, 'tag', tag_init));
        return; 
    end
    
    for i = 1:length(rowfields)
        for j = 1:length(colfields)

            row_fname = rowfields{i};
            col_fname = colfields{j};

            subplot(nrows, ncols, cnt_splt); 
            hold on; cnt_splt = cnt_splt + 1; 
            set(gca, 'Tag', row_fname); 

            plot(vis_results(j,1).res.(row_fname), '-', 'color', line_color, 'tag', tag_init);
            plot(vis_results(j,2).res.(row_fname), ':', 'color', line_color, 'tag', tag_init);
                
            if i == 1
                title(regexprep(col_fname, '_', '-'));
            end
            if j == 1
                ylabel(row_fname);
            end

        end
    end
    
    ylim_axes = struct; 
    ylim_axes.a = rng_of_a;
    ylim_axes.b = rng_of_b;    
    ylim_axes.dYcomp = [0, 1]; 
    ylim_axes.dYincomp = [0, 1]; 
    
    cellfun(@(tagname) linkaxes(findall(gcf, 'type', 'axes', 'tag', tagname), 'xy'), rowfields);
    cellfun(@(tagname) ylim(findall(gcf, 'type', 'axes', 'tag', tagname), ylim_axes.(tagname)), rowfields);
    if length(findobj(gca,'type','line')) == 2
        despline('all');
    end

    end