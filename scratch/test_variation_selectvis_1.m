clc; clear; close all; 
run start_up; 
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

train_results = train_results(select_cond_ind,:);
num_train_subconds = length(select_cond_ind); 
num_train_p_inc = length(train_p_inc_vec);

ind_t_select = ceil( size(train_results{1,1}.a.mean,2) * [1/2,1] ); 

bar_alpha_W_select = cellfun(@(x) x.bar_alpha_W.mean(:,ind_t_select), train_results, 'uni', 0); 
bar_alpha_W_select = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), bar_alpha_W_select, 'uni', 0));

cmap = [0.2,0.2,0.2; 0.2,0.3,0.8; 0.9,0.2,0.1; 0.2,0.5,0.1];

lgnd_names = {'no IP', 'only bias', 'only gain', 'both IP'};
figure('defaultaxesfontsize', 25, 'position', [0,0,0.8,1]);
subplot(311); hold on; 
hold on;colororder(cmap);
b2 = bar(train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,2)), 'Facealpha', 0.2);
arrayfun(@(h) set(h, 'edgecolor', h.FaceColor, 'linewidth', 2, 'facecolor', [0.8,0.8,0.8]), b2);

b1 = bar(train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,1)), 0.6,'FaceAlpha', 0.9,  'linestyle', 'none');
title('example summary plot (inside bars $\sim 0.5T$, outside $\sim T$)');
legend(b1, lgnd_names(:), 'NumColumns', 1, 'Location', 'west', 'fontsize', 20);
ylabel('$\overline{\alpha_W}$');
subplot(312); hold on;
hold on;colororder(cmap);
bar( train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,2)), 'FaceAlpha', 0.8, 'linestyle', 'none');
b1 = bar( train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,1)), 0.3, 'FaceAlpha', 0.9); 

arrayfun(@(h) set(h, 'facecolor', [1,1,1], 'edgecolor', [1,1,1], 'linestyle', 'none'), b1);

subplot(313); hold on; 
hold on;colororder(cmap); 
bar( train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,2)), 'FaceAlpha', 0.5, 'linestyle', 'none');
bar( train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,1)),0.4,'FaceAlpha', 0.8,  'linestyle', 'none'); 
xlabel('$p^{\mathrm{train}}_{\mathrm{inc}}$');
despline('all', [0.1,1.5]);
arrayfun(@(x) set(x.BaseLine,'linestyle', 'none'), findall(gcf, 'type', 'bar'))

%%
ind_t_select = ceil( size(train_results{1,1}.a.mean,2) * [1/2,1] ); 

bar_alpha_W_select = cellfun(@(x) x.bar_alpha_W.mean(:,ind_t_select), train_results, 'uni', 0); 
bar_alpha_W_select = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), bar_alpha_W_select, 'uni', 0));

dY_objvsnoise_test_select = cellfun(@(x) x.dY_objvsnoise_test.mean([1,3],ind_t_select), train_results, 'uni', 0);
dY_objvsnoise_test_select = cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), dY_objvsnoise_test_select, 'uni', 0));

cmap = [0.5,0.5,0.5; 0.2,0.3,0.8; 0.9,0.2,0.1; 0.8,0.2,0.7];
linestyles = {'--o', '-s'}; 
commonstyles = {'linewidth', 2, 'markersize', 20};

lgnd_names = {'no IP', 'only bias', 'only gain', 'both IP'};
figure; 
subplot(311); hold on; 
hold on;colororder(cmap);
b2 = bar(train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,2)), 'Facealpha', 0.2);
arrayfun(@(h) set(h, 'edgecolor', h.FaceColor, 'linewidth', 2, 'facecolor', [0.8,0.8,0.8]), b2);

bar(train_p_inc_vec, squeeze(bar_alpha_W_select(:,:,:,1)), 0.6,'FaceAlpha', 0.9,  'linestyle', 'none');
legend(b2, lgnd_names(:), 'NumColumns', 4, 'Location', 'northoutside');

subplot(312); hold on;
hold on;colororder(cmap);
b2 = bar( train_p_inc_vec, squeeze(dY_objvsnoise_test_select(:,:,1,2)), 'FaceAlpha', 0.8, 'linestyle', 'none');
b1 = bar( train_p_inc_vec, squeeze(dY_objvsnoise_test_select(:,:,1,1)), 0.3, 'FaceAlpha', 0.9); 

arrayfun(@(h) set(h, 'facecolor', [1,1,1], 'edgecolor', [1,1,1], 'linestyle', 'none'), b1);


subplot(313); hold on; 
hold on;colororder(cmap); 
bar( train_p_inc_vec, squeeze(dY_objvsnoise_test_select(:,:,2,2)), 'FaceAlpha', 0.5, 'linestyle', 'none');
bar( train_p_inc_vec, squeeze(dY_objvsnoise_test_select(:,:,2,1)),0.4,'FaceAlpha', 0.8,  'linestyle', 'none'); 

despline('all');
arrayfun(@(x) set(x.BaseLine,'linestyle', 'none'), findall(gcf, 'type', 'bar'))
%%
a_select = cellfun(@(x) x.a.mean, train_results, 'uni', 0); 
a_select = permute(cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), a_select, 'uni', 0)), [2,3,1,4]);
a_select = reshape(a_select, [size(a_select,[1,2]), prod(size(a_select,[3,4]))]);

b_select = cellfun(@(x) x.b.mean, train_results, 'uni', 0); 
b_select = permute(cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), b_select, 'uni', 0)), [2,3,1,4]);
b_select = reshape(b_select, [size(b_select,[1,2]), prod(size(b_select,[3,4]))]);

bar_alpha_W_select = cellfun(@(x) x.bar_alpha_W.mean, train_results, 'uni', 0); 
bar_alpha_W_select = permute(cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), bar_alpha_W_select, 'uni', 0)), [2,3,1,4]);
bar_alpha_W_select = reshape(bar_alpha_W_select, [size(bar_alpha_W_select,[1,2]), prod(size(bar_alpha_W_select,[3,4]))]);

dY_objvsnoise_test_select = cellfun(@(x) x.dY_objvsnoise_test.mean, train_results, 'uni', 0); 
dY_objvsnoise_test_select = permute(cell2mat(cellfun(@(x) reshape(x, [1,1,size(x)]), dY_objvsnoise_test_select, 'uni', 0)), [2,3,1,4]);
dY_objvsnoise_test_select = reshape(dY_objvsnoise_test_select, [size(dY_objvsnoise_test_select,[1,2]), prod(size(dY_objvsnoise_test_select,[3,4]))]);

rhos = arrayfun(@(i) corr(squeeze(cat(2,a_select(i,:,:),b_select(i,:,:)))', ...
    squeeze(cat(2,bar_alpha_W_select(i,:,:),dY_objvsnoise_test_select(i,:,:)))'), 1:6, 'uni', 0);
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

% fig_name = sprintf('%s_properties', file_prefix); 
% export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
% close; 


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

% fig_name = sprintf('%s_perfproxy', file_prefix); 
% export_fig(fullfile(fig_path, fig_name), '-r300', '-p0.02'); 
% close; 
