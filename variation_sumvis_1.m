clc; clear; close all; 

file_prefix = 'full_var-IP-onoff_var-norminput-onoff_IPfasterthanSP'; 
load(fullfile('data', file_prefix), ...
    'def_opts', 'train_p_inc_vec', 'train_results');


fig_path = 'figures/sum_1'; 
if ~exist(fig_path, 'dir') 
    mkdir(fig_path);
end

%% Select only `train_results` with normalized inputs
select_cond = cellfun(@(x) x.opts.train_norm_input, train_results);
select_cond_ind = find(sum(select_cond,2) == length(train_p_inc_vec)); 

num_train_subconds = length(select_cond_ind); 
num_train_p_inc = length(train_p_inc_vec);

latex_ips = {['no IP' newline '($\eta_{a,b} = 0$)'], ['only bias plast.' newline '(only $\eta_b > 0$)'];
            ['only gain plast.' newline '(only $\eta_a > 0$)'], ['full IP' newline '($\eta_{a,b} > 0$)']};
        
%% Plot Properties
plt_fields = {'a', 'b', 'bar_alpha_W'}; 
latex_fields = {'a \sim \mathrm{gain}','b \sim \mathrm{bias}','\overline{\alpha_W}'};
nrows = length(plt_fields); 
ncols = num_train_subconds; 

step_vec = 1:def_opts.num_train;
step_vec = step_vec(1:def_opts.subsampled:end);

train_transition = ceil(def_opts.num_train * def_opts.p_train_complete); 

cmap = parula(num_train_p_inc)*0.85; 

figure('WindowStyle', 'normal', ...
    'defaultaxeslabelfontsize', 1.4, ...
    'defaultaxestitlefontsize', 1.4, ...
    'Position', [0,0,1,0.7]);

cnt_splt = 1; 
for i = 1:nrows
    for j = 1:ncols
        subplot(nrows, ncols, cnt_splt); hold on; 
        set(gca, 'tag', plt_fields{i}); 
        cnt_splt = cnt_splt + 1; 
        
        plt_field = plt_fields{i}; 
        res_row_ind = select_cond_ind(j);
        
        dat2plt = vertcat(vertcat(train_results{res_row_ind,:}).(plt_field));
        arrayfun(@(x) plot(step_vec, dat2plt(x).mean, ...
            'linewidth', 2, 'color', cmap(x,:)), 1:num_train_p_inc); 
        arrayfun(@(x) fill([step_vec, fliplr(step_vec)], ...
            [(dat2plt(x).mean + dat2plt(x).sem), fliplr(dat2plt(x).mean - dat2plt(x).sem)], ...
            cmap(x,:), 'LineStyle', 'none', 'FaceColor', cmap(x,:), 'FaceAlpha', 0.3), 1:num_train_p_inc);
        
        if i == 1
            train_subcond_opts = train_results{res_row_ind,1}.opts;

            subplot_description = latex_ips{...
                (train_subcond_opts.eta_ip_a > 0) + 1, ...
                (train_subcond_opts.eta_ip_b > 0) + 1};
            title(subplot_description);
        end
        
        if j == 1
            ylabel(sprintf('$%s$', latex_fields{i}));            
        end
        
        if i == nrows
            xlabel('\# step');
        else
            set(gca, 'xcolor', 'none');
        end
        
        if i == 1 && j == ncols
            cbar = colorbar;
            colormap(cmap);
            caxis(train_p_inc_vec([1,end]));
            title(cbar, '$p_{\mathrm{inc}}^{\mathrm{train}}$', 'interpreter', 'latex');
            cbar.Position = cbar.Position .* [1,1,1,0.7] + [0.08,0.02,0,0];
            cbar.FontSize = 25;
        end        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
ylim(findall(gcf, 'type', 'axes', 'tag', 'a'), [5,20]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'b'), [-5,-3.5]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'bar_alpha_W'), [-0.1, 0.9]);

cellfun(@(s) linkaxes(findall(gcf, 'type', 'axes', 'tag', s), 'y'), plt_fields);
linkaxes(findall(gcf, 'type', 'axes'), 'x');
despline('all');
xlim([0,1]*def_opts.num_train);

fig_name = sprintf('%s_properties', file_prefix); 
export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
close; 


%% Plot Performance Proxy

latex_ips = {'no IP', 'only bias plast.';
            'only gain plast.', 'full IP'};
 
select_test_p_inc_ind = [1,3]; 
latex_fields = arrayfun(@(x) sprintf('p_{\\mathrm{inc}}^{\\mathrm{test}} = %g', x),...
    def_opts.test_p_incs(select_test_p_inc_ind), 'uni', 0); 
nrows = length(latex_fields); 
ncols = num_train_subconds; 

figure('WindowStyle', 'normal', ...
    'units', 'normalized', 'Position', [0,0,1,0.7], ...
    'defaultaxeslabelfontsize', 1.4, ...
    'defaultaxestitlefontsize', 1.4);

cnt_splt = 1; 
for i = 1:nrows
    for j = 1:ncols
        subplot(nrows, ncols, cnt_splt); hold on; 
        cnt_splt = cnt_splt + 1; 
        
        res_row_ind = select_cond_ind(j);
        stat_row_ind = select_test_p_inc_ind(i);
        
        dat2plt = vertcat(vertcat(train_results{res_row_ind,:}).dY_objvsnoise_test);
        
        arrayfun(@(x) plot(step_vec, dat2plt(x).mean(stat_row_ind,:), ...
            'linewidth', 2, 'color', cmap(x,:)), 1:num_train_p_inc); 
        arrayfun(@(x) fill([step_vec, fliplr(step_vec)], ...
            [(dat2plt(x).mean(stat_row_ind,:) + dat2plt(x).sem(stat_row_ind,:)), ...
            fliplr(dat2plt(x).mean(stat_row_ind,:) - dat2plt(x).sem(stat_row_ind,:))], ...
            cmap(x,:), 'LineStyle', 'none', 'FaceColor', cmap(x,:), 'FaceAlpha', 0.3), 1:num_train_p_inc);
        
        
        if i == 1
    
            train_subcond_opts = train_results{res_row_ind,1}.opts;
         
            subplot_description = latex_ips{...
                (train_subcond_opts.eta_ip_a > 0) + 1, ...
                (train_subcond_opts.eta_ip_b > 0) + 1};
            title(subplot_description);
            
            title(subplot_description);
        end
        
        if j == 1
            % ylabel(['$\langle Y_{\mathrm{selective}}\rangle - \langle Y_{\mathrm{noise}}\rangle$' newline ...
            %    sprintf('($%s$)', latex_fields{i})]);           
            % ylabel(['$\Delta Y^{\mathrm{(test)}}_{\mathrm{selec - noise}} \vert ' latex_fields{i} '$']);
            ylabel(['$\Delta Y^{\mathrm{(test)}} \vert ' latex_fields{i} '$']);
        end
        
        if i == nrows
            xlabel('\# step');
        else
            set(gca, 'xcolor', 'none');
        end
        
        if i == nrows && j == ncols
            cbar = colorbar;
            colormap(cmap);
            caxis(train_p_inc_vec([1,end]));
            title(cbar, '$p_{\mathrm{inc}}^{\mathrm{train}}$', 'interpreter', 'latex');
            cbar.Position = cbar.Position .* [1,1,1,0.5] + [0.06,0.02,0,0];
            cbar.FontSize = 25;
        end
        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
linkaxes(findall(gcf, 'type', 'axes'), 'xy');

xlim([0,1]*def_opts.num_train);
ylim([0,1]);
despline('all',[1e-1,2e-1]);

fig_name = sprintf('%s_perfproxy', file_prefix); 
export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
close; 


%% Select only `train_results` with normalized inputs
barplot_style = 2;
select_test_p_inc_ind = [1,3]; 

select_cond = cellfun(@(x) x.opts.train_norm_input, train_results);
select_cond_ind = find(sum(select_cond,2) == length(train_p_inc_vec)); 

train_results = train_results(select_cond_ind,:);
num_train_subconds = length(select_cond_ind); 
num_train_p_inc = length(train_p_inc_vec);

ind_t_select = ceil( size(train_results{1,1}.a.mean,2) * [1/2,1] ); 

bar_alpha_W_select = cellfun(@(x) x.bar_alpha_W.mean(:,ind_t_select), train_results, 'uni', 0); 
bar_alpha_W_select = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), bar_alpha_W_select, 'uni', 0));

dY_objvsnoise_test_select = cellfun(@(x) x.dY_objvsnoise_test.mean(select_test_p_inc_ind,ind_t_select), train_results, 'uni', 0);
dY_objvsnoise_test_select = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), dY_objvsnoise_test_select, 'uni', 0));

data_barplot = struct;
data_barplot.alpha.ylabel = '$\overline{\alpha_W}$';
data_barplot.alpha.title = 'summary of selectivity formation (inside bars $\sim 0.5T$, outside $\sim T$)';
data_barplot.alpha.data.during_training = squeeze(bar_alpha_W_select(:,:,:,1));
data_barplot.alpha.data.end_training = squeeze(bar_alpha_W_select(:,:,:,2));

dYtest_latex_string = '\Delta Y^{\mathrm{(test)}} \vert p_{\mathrm{inc}}^{\mathrm{test}}';
data_barplot.dYtestcomp.ylabel = sprintf('$%s=%g$', dYtest_latex_string, def_opts.test_p_incs(select_test_p_inc_ind(1)));
data_barplot.dYtestcomp.title = 'summary of performance proxy, when tested with \textbf{complete} input';
data_barplot.dYtestcomp.data.during_training = squeeze(dY_objvsnoise_test_select(:,:,1,1));
data_barplot.dYtestcomp.data.end_training = squeeze(dY_objvsnoise_test_select(:,:,1,2));

data_barplot.dYtestincomp.ylabel = sprintf('$%s=%g$', dYtest_latex_string, def_opts.test_p_incs(select_test_p_inc_ind(2)));
data_barplot.dYtestincomp.title = 'summary of performance proxy, when tested with \textbf{incomplete} input';
data_barplot.dYtestincomp.data.during_training = squeeze(dY_objvsnoise_test_select(:,:,2,1));
data_barplot.dYtestincomp.data.end_training = squeeze(dY_objvsnoise_test_select(:,:,2,2));

cmap_barplot = [0.2,0.2,0.2; 0.2,0.3,0.8; 0.9,0.2,0.1; 0.2,0.5,0.1];

lgnd_names = {'no IP', 'only bias plast.', 'only gain plast.', 'full IP'};
barplot_fields = fieldnames(data_barplot); 

figure('defaultaxesfontsize', 18, ...
    'WindowStyle', 'normal', ...
    'units', 'normalized', 'Position', [0,0,0.8,1], ...
    'defaultaxeslabelfontsize', 1.2, ...
    'defaultaxestitlefontsize', 1.2);

for i = 1:length(barplot_fields)
    subplot(length(barplot_fields),1,i); hold on;
    colororder(cmap_barplot); 
    
    field_ith = barplot_fields{i};
    data_ith = data_barplot.(field_ith); 
    
    switch barplot_style
        case 1
            bar_endtraining = bar(train_p_inc_vec, data_ith.data.end_training, 'Facealpha', 0.2);
            arrayfun(@(h) set(h, 'edgecolor', h.FaceColor, 'linewidth', 2, 'facecolor', [0.8,0.8,0.8]), bar_endtraining);
            
            bar_duringtraining = bar(train_p_inc_vec, data_ith.data.during_training, 0.6,'FaceAlpha', 0.9,  'linestyle', 'none');
            lgnd_bar_handle = bar_duringtraining; 
        case 2
            bar_endtraining = bar(train_p_inc_vec, data_ith.data.end_training, 'FaceAlpha', 0.8, 'linestyle', 'none');
            bar_duringtraining = bar(train_p_inc_vec, data_ith.data.during_training, 0.3, 'FaceAlpha', 0.9);
            
            arrayfun(@(h) set(h, 'facecolor', [1,1,1], 'edgecolor', [1,1,1], 'linestyle', 'none'), bar_duringtraining);
            lgnd_bar_handle = bar_endtraining;
        case 3            
            bar(train_p_inc_vec, data_ith.data.end_training, 'FaceAlpha', 0.5, 'linestyle', 'none');
            bar(train_p_inc_vec, data_ith.data.during_training, 0.4,'FaceAlpha', 0.8,  'linestyle', 'none');
            lgnd_bar_handle = bar_endtraining;
        otherwise 
            error('not allowed');
    end
    ylabel(data_ith.ylabel); 
    title(data_ith.title); 
    
    ax = gca; 
    ax.Position = ax.Position + [0,-0.02,0,0];
    
    if i < length(barplot_fields)
        ax.XColor = 'none';
    end
    
    if i == 1
        lgnd_obj = legend(lgnd_bar_handle, lgnd_names(:), 'NumColumns', length(lgnd_names), 'Location', 'north', 'fontsize', 18, 'box', 'on');
        lgnd_obj.Position = lgnd_obj.Position + [0,0.1,0,0];
    end
        
end

xlabel('input training incompleteness ($p^{\mathrm{train}}_{\mathrm{inc}}$)');
despline('all', [0.1,1.5]);
arrayfun(@(x) set(x.BaseLine,'linestyle', 'none'), findall(gcf, 'type', 'bar'))

fig_name = sprintf('%s_summary-select', file_prefix); 
export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
close; 


