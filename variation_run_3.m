clc; clear; close all; 
run start_up; 

file_prefix = 'full_var-IP-onoff_IP-test3'; 
max_eta_ip = 1e-2;

%% Default opts 
def_opts = struct; 
def_opts.N = 10;

def_opts.mu = 1/(2*def_opts.N);
def_opts.selrow = 8; 

def_opts.num_train = 10e4; 
def_opts.p_train_complete = 0;

def_opts.eta_sp_ltp = 1e-3; 
def_opts.kappa_sp = 1/20; % kappa_sp=eta_sp_ltd/eta_sp_ltp = 0, 1/20, 1/15, 1/10

def_opts.weight_properties.norm_order = 2; % 1 or 2
def_opts.weight_properties.onlyexcitatory = true; % true or false 

def_opts.subsampled = 1000;

def_opts.num_test_per_pinc = 1000;
def_opts.test_p_sel = 0.5; % 50-50 of selective input 
def_opts.test_sigma_noise = 0.1; % small gaussian noise 
def_opts.test_p_incs = [0, 0.1, 0.3, 0.5]; 
def_opts.test_shuffle = false; % not necessary 

def_opts.k_r = 1; 
def_opts.norm_input_opt = 1; 
def_opts.a_init = 4; 
def_opts.b_init = 0; 


%% Variations 
train_p_inc_vec = [0,0.1,0.3,0.5];
num_rand = 10; 

% for eta_ip_a, eta_ip_b
eta_ip_vec = [0,max_eta_ip];

% create varitions
[eta_a_s, eta_b_s] = meshgrid(eta_ip_vec, eta_ip_vec); 

training_subconditions = struct(...
    'eta_ip_a', num2cell(eta_a_s(:)), ...
    'eta_ip_b', num2cell(eta_b_s(:)));
vec_fields_subcond = fieldnames(training_subconditions); 

num_train_p_inc_vec = length(train_p_inc_vec);
num_train_subconds = length(training_subconditions);

%% Train 

train_results = cell(num_train_subconds, num_train_p_inc_vec); 
vec_fields_tosave = {'a','b','bar_alpha_W', 'vec_train_completeness', 'dY_objvsnoise_test'};

for i = 1:num_train_subconds
    train_cond_struct = training_subconditions(i);
    
    tic
    tmp_res = cell(1, num_train_p_inc_vec); 
    parfor j = 1:num_train_p_inc_vec
        train_p_inc = train_p_inc_vec(j); 
        
        opts = def_opts;
        opts.train_p_inc = train_p_inc;
        for k = 1:length(vec_fields_subcond)
            opts.(vec_fields_subcond{k}) = train_cond_struct.(vec_fields_subcond{k}); 
        end
        
        res_j = arrayfun(@(~) train_with_incomplete(opts), 1:num_rand, 'uni', 1);
        tmp_res_j = struct; 
        for k = 1:length(vec_fields_tosave)
            cat_vec = cat(1,res_j.(vec_fields_tosave{k}));
            tmp_res_j.(vec_fields_tosave{k}).mean = squeeze(mean(cat_vec));
            tmp_res_j.(vec_fields_tosave{k}).sem = squeeze(std(cat_vec) / sqrt(num_rand));
        end   
        
        tmp_res_j.opts = opts;
        tmp_res{j} = tmp_res_j;
    end
    toc
    
    train_results(i,:) = tmp_res;
end

save(fullfile('data', file_prefix), ...
    'def_opts', 'train_p_inc_vec', 'training_subconditions', 'train_results');

%% Plot Properties
plt_fields = {'a', 'b', 'bar_alpha_W', 'vec_train_completeness'}; 
latex_fields = {'a','b','\overline{\alpha_W}', '\%c'};
nrows = length(plt_fields); 
ncols = num_train_subconds; 

step_vec = 1:def_opts.num_train;
step_vec = step_vec(1:def_opts.subsampled:end);

train_transition = ceil(def_opts.num_train * def_opts.p_train_complete); 

cmap = parula(num_train_p_inc_vec)*0.85; 

var_description = sprintf('[$\\eta_a$ on (%g)/off, $\\eta_b$ on (%g)/off, norm input]', ...
    max(eta_ip_vec), max(eta_ip_vec));
fig_description = sprintf('$\\eta_{sp}^{(LTP)} = %g, \\kappa_{sp} = %g, N = %d$ %s', ...
    def_opts.eta_sp_ltp, def_opts.kappa_sp, def_opts.N, var_description); 

figure; 

annotation('textbox', 'units', 'normalized', 'position', [0,0.92,1,0.08], ...
    'string', fig_description, 'linestyle', 'none', ...
    'fontsize', 15, 'interpreter', 'latex', ...
    'horizontalalignment', 'center', 'verticalalignment', 'top'); 

cnt_splt = 1; 
for i = 1:nrows
    for j = 1:ncols
        subplot(nrows, ncols, cnt_splt); hold on; 
        set(gca, 'tag', plt_fields{i}); 
        cnt_splt = cnt_splt + 1; 
        
        plt_field = plt_fields{i}; 
        
        dat2plt = vertcat(vertcat(train_results{j,:}).(plt_field));
        arrayfun(@(x) plot(step_vec, dat2plt(x).mean, ...
            'linewidth', 2, 'color', cmap(x,:)), 1:num_train_p_inc_vec); 
        arrayfun(@(x) fill([step_vec, fliplr(step_vec)], ...
            [(dat2plt(x).mean + dat2plt(x).sem), fliplr(dat2plt(x).mean - dat2plt(x).sem)], ...
            cmap(x,:), 'LineStyle', 'none', 'FaceColor', cmap(x,:), 'FaceAlpha', 0.3), 1:num_train_p_inc_vec);
        
        if i == 1
            train_subcond_opts = train_results{j,1}.opts;
            subplot_description = sprintf(...
                '[%d, %d]', ...
                train_subcond_opts.eta_ip_a > 0, ...
                train_subcond_opts.eta_ip_b > 0);

            title(subplot_description);
        end
        
        if j == 1
            ylabel(sprintf('$%s$', latex_fields{i}));            
        end
        
        if i == nrows
            xlabel('\# step');
        end
        
        if i == nrows && j == ncols
            cbar = colorbar;
            colormap(cmap);
            caxis(train_p_inc_vec([1,end]));
            title(cbar, '$p_{inc}^{train}$', 'interpreter', 'latex');
            cbar.Position = cbar.Position .* [1,1,1,0.5] + [0.08,0.02,0,0];
        end
        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
% ylim(findall(gcf, 'type', 'axes', 'tag', 'a'), [3,25]);
% ylim(findall(gcf, 'type', 'axes', 'tag', 'b'), [-6,-3]);
cellfun(@(s) linkaxes(findall(gcf, 'type', 'axes', 'tag', s), 'y'), plt_fields);
linkaxes(findall(gcf, 'type', 'axes'), 'x');
despline('all');
xlim([0.01,1]*def_opts.num_train);

% fig_name = sprintf('%s_properties', file_prefix);
% savefig(fullfile('figures', fig_name));
%  
% export_fig(fullfile('figures', fig_name), '-r300', '-p0.02'); 
% close; 


%% Plot Performance Proxy

latex_fields = arrayfun(@(x) sprintf('p_{inc}^{test} = %g', x),...
    def_opts.test_p_incs, 'uni', 0); 
nrows = length(latex_fields); 
ncols = num_train_subconds; 

figure; 

annotation('textbox', 'units', 'normalized', 'position', [0,0.92,1,0.08], ...
    'string', fig_description, 'linestyle', 'none', ...
    'fontsize', 20, 'interpreter', 'latex', ...
    'horizontalalignment', 'center', 'verticalalignment', 'top'); 

cnt_splt = 1; 
for i = 1:nrows
    for j = 1:ncols
        subplot(nrows, ncols, cnt_splt); hold on; 
        cnt_splt = cnt_splt + 1; 
        
        dat2plt = vertcat(vertcat(train_results{j,:}).dY_objvsnoise_test);
        arrayfun(@(x) plot(step_vec, dat2plt(x).mean(i,:), ...
            'linewidth', 2, 'color', cmap(x,:)), 1:num_train_p_inc_vec); 
        arrayfun(@(x) fill([step_vec, fliplr(step_vec)], ...
            [(dat2plt(x).mean(i,:) + dat2plt(x).sem(i,:)), fliplr(dat2plt(x).mean(i,:) - dat2plt(x).sem(i,:))], ...
            cmap(x,:), 'LineStyle', 'none', 'FaceColor', cmap(x,:), 'FaceAlpha', 0.3), 1:num_train_p_inc_vec);
        
        
        if i == 1
    
            train_subcond_opts = train_results{j,1}.opts;
            subplot_description = sprintf(...
                '[%d, %d]', ...
                train_subcond_opts.eta_ip_a > 0, ...
                train_subcond_opts.eta_ip_b > 0);
            
            title(subplot_description);
        end
        
        if j == 1
            ylabel(sprintf('$%s$', latex_fields{i}));            
        end
        
        if i == nrows
            xlabel('\# step');
        end
        
        if i == nrows && j == ncols
            cbar = colorbar;
            colormap(cmap);
            caxis(train_p_inc_vec([1,end]));
            title(cbar, '$p_{inc}^{train}$', 'interpreter', 'latex');
            cbar.Position = cbar.Position .* [1,1,1,0.5] + [0.08,0.02,0,0];
        end
        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
linkaxes(findall(gcf, 'type', 'axes'), 'xy');
despline('all');
xlim([0.01,1]*def_opts.num_train);
ylim([0,1]);

% fig_name = sprintf('%s_perfproxy', file_prefix);
% savefig(fullfile('figures', fig_name));
%  
% export_fig(fullfile('figures', fig_name), '-r300', '-p0.02'); 
% close; 

