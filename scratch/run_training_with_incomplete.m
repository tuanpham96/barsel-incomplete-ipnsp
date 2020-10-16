clc; clear; close all; 
run start_up; 

file_prefix = 'run_incomplete_scan_eta_ip_a_and_p_inc_fastsp'; 

%% Default opts 
def_opts = struct; 
def_opts.N = 10;
def_opts.a_init = 5;
def_opts.b_init = -5;
def_opts.mu = 1/(2*def_opts.N);
def_opts.selrow = 8; 

def_opts.num_train = 100e4; 
def_opts.p_train_complete = 1/10;

def_opts.eta_sp = 1e-2; 
def_opts.eta_ip_b = 1e-2;

def_opts.subsampled = 100; 


%% Variations 
p_inc_vec = 0.1:0.2:0.9;
eta_ip_a_vec = [0,1e-3,1e-2,1e-1];

num_rand = 5; 

num_p_inc_vec = length(p_inc_vec);
num_eta_ip_a_vec = length(eta_ip_a_vec); 

train_results = cell(num_eta_ip_a_vec, num_p_inc_vec); 
vec_fields = {'a','b','alphaW'};

for i = 1:num_eta_ip_a_vec
    eta_ip_a = eta_ip_a_vec(i);
    
    tic
    tmp_res = cell(1, num_p_inc_vec); 
    parfor j = 1:num_p_inc_vec
        p_inc = p_inc_vec(j); 
        
        opts = def_opts;
        opts.p_inc = p_inc;
        opts.eta_ip_a = eta_ip_a;
        
        res_j = arrayfun(@(~) train_with_incomplete(opts), 1:num_rand, 'uni', 1);
        tmp_res_j = struct; 
        for k = 1:length(vec_fields)
            cat_vec = vertcat(res_j.(vec_fields{k}));
            tmp_res_j.(vec_fields{k}).mean = mean(cat_vec);
            tmp_res_j.(vec_fields{k}).sem = std(cat_vec) / sqrt(num_rand);
        end
        
        
        tmp_res{j} = tmp_res_j;
    end
    toc
    
    train_results(i,:) = tmp_res;
end

save(fullfile('data', file_prefix), ...
    'def_opts', 'p_inc_vec', 'eta_ip_a_vec', 'train_results');

%% Plot
plt_fields = {'a', 'b', 'alphaW'}; 
latex_fields = {'a','b','\alpha_{norm-min-max(W)}'};
nrows = length(plt_fields); 
ncols = num_eta_ip_a_vec; 

step_vec = 1:def_opts.num_train;
step_vec = step_vec(1:def_opts.subsampled:end);

train_transition = ceil(def_opts.num_train * def_opts.p_train_complete); 

cmap = parula(num_p_inc_vec)*0.85; 
fig_description = sprintf('$\\eta_{ip,b} = %g, \\eta_{sp} = %g, N = %d$', ...
    def_opts.eta_ip_b, def_opts.eta_sp, def_opts.N); 

figure; 

annotation('textbox', 'units', 'normalized', 'position', [0,0.92,1,0.08], ...
    'string', fig_description, 'linestyle', 'none', ...
    'fontsize', 20, 'interpreter', 'latex', ...
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
            'linewidth', 2, 'color', cmap(x,:)), 1:num_p_inc_vec); 
        arrayfun(@(x) fill([step_vec, fliplr(step_vec)], ...
            [(dat2plt(x).mean + dat2plt(x).sem), fliplr(dat2plt(x).mean - dat2plt(x).sem)], ...
            cmap(x,:), 'LineStyle', 'none', 'FaceColor', cmap(x,:), 'FaceAlpha', 0.3), 1:num_p_inc_vec);
        
        % title(sprintf('$%s, \\eta_{ip,a} = %g$', latex_fields{i}, eta_ip_a_vec(j)));
        
        if i == 1
            title(sprintf('$\\eta_{ip,a} = %g$', eta_ip_a_vec(j)));
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
            caxis(p_inc_vec([1,end]));
            title(cbar, 'p_{inc}');
            cbar.Position = cbar.Position .* [1,1,1,0.5] + [0.08,0.02,0,0];
        end
        
       
        xline(train_transition, ':k', 'linewidth', 1); 
        
    end
end
ylim(findall(gcf, 'type', 'axes', 'tag', 'a'), [3,25]);
ylim(findall(gcf, 'type', 'axes', 'tag', 'b'), [-6,-3]);
cellfun(@(s) linkaxes(findall(gcf, 'type', 'axes', 'tag', s), 'y'), plt_fields);
linkaxes(findall(gcf, 'type', 'axes'), 'x'),
despline('all');
xlim([0.01,1]*def_opts.num_train);

export_fig(fullfile('figures', file_prefix), '-r300', '-p0.02'); 
close; 