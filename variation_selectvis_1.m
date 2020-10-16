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

latex_ips = {'no IP', 'only $\eta_b > 0$';
            'only $\eta_a > 0$', 'IP with both $a,b$'};
        
%% Plot Properties
plt_fields = {'a', 'b', 'bar_alpha_W', 'vec_train_completeness'}; 
latex_fields = {'a','b','\overline{\alpha_W}', '\%c'};
nrows = length(plt_fields); 
ncols = num_train_subconds; 

step_vec = 1:def_opts.num_train;
step_vec = step_vec(1:def_opts.subsampled:end);

train_transition = ceil(def_opts.num_train * def_opts.p_train_complete); 

cmap = parula(num_train_p_inc)*0.85; 

figure('WindowStyle', 'normal', ...
    'defaultaxeslabelfontsize', 1.6, ...
    'defaultaxestitlefontsize', 1.6);

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
        
        if i == nrows
            ax = gca; 
            ax.Position = ax.Position .* [1,1,1,0.5] + [0,0.1,0,0];
        end
        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
ylim(findall(gcf, 'type', 'axes', 'tag', 'a'), [0,20]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'b'), [-5,-3.5]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'bar_alpha_W'), [-0.1, 0.8]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'vec_train_completeness'), [0.4, 1.1]);

cellfun(@(s) linkaxes(findall(gcf, 'type', 'axes', 'tag', s), 'y'), plt_fields);
linkaxes(findall(gcf, 'type', 'axes'), 'x');
despline('all');
xlim([0,1]*def_opts.num_train);

fig_name = sprintf('%s_properties', file_prefix); 
export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
close; 


%% Plot Performance Proxy
select_test_p_inc_ind = [1,3]; 
latex_fields = arrayfun(@(x) sprintf('p_{\\mathrm{inc}}^{\\mathrm{test}} = %g', x),...
    def_opts.test_p_incs(select_test_p_inc_ind), 'uni', 0); 
nrows = length(latex_fields); 
ncols = num_train_subconds; 

figure('WindowStyle', 'normal', ...
    'units', 'normalized', 'Position', [0,0,1,0.8], ...
    'defaultaxeslabelfontsize', 1.6, ...
    'defaultaxestitlefontsize', 1.6);

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
            ylabel(['$\langle Y_{\mathrm{selective}}\rangle - \langle Y_{\mathrm{noise}}\rangle$' newline ...
                sprintf('($%s$)', latex_fields{i})]);            
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
