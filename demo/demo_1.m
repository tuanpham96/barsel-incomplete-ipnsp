clc; clear; close all; 
run start_up.m; 

fig_path = 'figures/demo'; 
if ~exist(fig_path, 'dir') 
    mkdir(fig_path);
end

%% Demo for `p_inc_train` (or `train_p_inc` in `opts`)

rng(957);
N = 10;
p_bar = 1/(2*N);
norm_input_opt = 1; 
num_input = 100;
num_toplot = 15;
p_inc_train_demo = [0, 0.1, 0.3, 0.5]; 
fig_prefix = 'demo-p_inc_train'; 

for i = 1:length(p_inc_train_demo)
    p_inc_train = p_inc_train_demo(i); 
    Us = arrayfun(@(~) reshape(...
        generate_bar_input_incomplete(N, p_bar, p_inc_train, norm_input_opt), ...
        [N,N]), 1:num_input, 'uni', 0);
    
    % filter out most of the zeros for demo purposes
    Us_nonzero = Us(cellfun(@(x) sum(x(:)) > 0, Us));
    Us_zeros = Us(cellfun(@(x) sum(x(:)) == 0, Us));
    Us = [Us_zeros(2) Us_nonzero];
    Us = Us(1:num_toplot);
    Us = Us(randperm(num_toplot));
    
    figure('WindowStyle', 'normal');
    
    fig_description = sprintf('$p_{\\mathrm{inc}}^{\\mathrm{train}} = %g$', p_inc_train);
    annotation('textbox', 'units', 'normalized', 'position', [0.1,0.92,1,0.08], ...
        'string', fig_description, 'linestyle', 'none', ...
        'fontsize', 50, 'interpreter', 'latex', ...
        'horizontalalignment', 'left', 'verticalalignment', 'top');
    colormap('gray');
    for j = 1:num_toplot
        subplot(3,5,j); hold on;
        image_with_strict_limits(Us{j}');
        daspect([1,1,1]); box on;
        set(gca, 'xtick', '', 'ytick', '', 'linewidth', 3);
        ax_pos = get(gca, 'position'); 
        ax_pos = ax_pos + [0,-0.07,0,0];
        set(gca, 'position', ax_pos);
    end
    
    fig_name = fullfile(fig_path, sprintf('%s_%02d', fig_prefix, i)); 
    export_fig(fig_name, '-r300', '-p0.02');
    close;
    
    
end

%% Demo for `p_inc_test` (or `test_p_inc` in `opts`)

rng(957);
N = 10;
selrow = 8; 
p_sel = 0.8; % for demo purpose
shuffle_opt = 0; 
norm_input_opt = 1; 
sigma_noise = 0.1; 
num_data_per_pinc = 5; 
p_inc_test_demo = [0, 0.1, 0.3, 0.5]; 
fig_prefix = 'demo-p_inc_test'; 

[Us, labels, p_inc_vec] = generate_test_input(num_data_per_pinc, N, p_sel, ...
    selrow, sigma_noise, p_inc_test_demo, norm_input_opt, shuffle_opt);

figure;

reordered_ind = reshape(1:size(Us,2), [num_data_per_pinc, length(p_inc_test_demo)])';

colormap('gray');
for i = 1:size(Us,2)
    subplot(4,5,reordered_ind(i)); hold on;
    image_with_strict_limits(reshape(Us(:,i),[N,N]));
    daspect([1,1,1]); box on;
    set(gca, 'xtick', '', 'ytick', '', 'linewidth', 3);
    if i <= length(p_inc_test_demo)
        panel_description = sprintf('$p_{\\mathrm{inc}}^{\\mathrm{test}} = %g$', p_inc_test_demo(i));
        ylabel(panel_description, 'fontsize', 30);
    end
end

fig_name = fullfile(fig_path, sprintf('%s', fig_prefix));
export_fig(fig_name, '-r300', '-p0.02');
close;

%% Demo plasticity 
rng(24829); 
opts = struct; 
opts.N = 10;
opts.a_init = 5;
opts.b_init = 5; 
opts.mu = 1/(2*opts.N);
opts.selrow = 8; 

opts.train_norm_input = true; 
opts.num_train = 30e4; 
opts.p_train_complete = 1/15;

opts.num_test_per_pinc = 1000;
opts.train_p_inc = 0.3;
opts.test_p_sel = 0.5;
opts.test_sigma_noise = 0.1; 
opts.test_norm_input = true; 
opts.test_p_incs = [0, 0.1, 0.3, 0.5]; 
opts.test_shuffle = false; % not necessary 

opts.eta_sp_ltp = 1e-3; 
opts.kappa_sp = 1/20; 
opts.weight_properties.norm_order = 2; 
opts.weight_properties.onlyexcitatory = true; 

opts.eta_ip_b = 1e-2;
opts.eta_ip_a = 1e-2; 

opts.subsampled = 500;

opts.points_to_record_W = [2e3,5e4,10e4,20e4,28e4]; 

tic
demo_res = train_with_incomplete(opts); 
toc

%% Some setup for graphics
select_steps = opts.points_to_record_W;
cmap = return_colorbrewer('Blues', length(select_steps)) * 0.9; 
ind_of_subsampled = arrayfun(@(x) find_nearest(demo_res.t, x, 'ind'), select_steps);
a_select = demo_res.a(ind_of_subsampled);
b_select = demo_res.b(ind_of_subsampled); 

%% Plot select `W`

fig_prefix = 'demo-training-select-W';
num_Ws_toplot = length(opts.points_to_record_W);
for i = 1:num_Ws_toplot    
    figure('WindowStyle', 'normal'); hold on; colormap('gray'); 
    image_with_strict_limits(demo_res.W{i}); 
    caxis([0, 0.35]);
    daspect([1,1,1]); box on;
    set(gca, 'xtick', '', 'ytick', '', 'linewidth', 3);
    
    fig_name = fullfile(fig_path, sprintf('%s_atstep%05d', fig_prefix, opts.points_to_record_W(i)));
    export_fig(fig_name, '-r300', '-p0.02');
    close;
    
end

%% Plot `a`, `b`, `bar_alpha_W`
fig_prefix = 'demo-training-a-b-baralpha';
train_transition = ceil(opts.num_train * opts.p_train_complete); 

lgdn_objs = gobjects(1,3); 

figure('WindowStyle', 'normal', 'DefaultLineLineWidth', 4, ...
    'Units', 'normalized', 'Position', [0,0.1,1,0.8]); 

subplot(311); hold on;
plot(demo_res.t, demo_res.vec_train_completeness, '-k'); 
set(gca, 'xcolor', 'none');
ylabel(['train-input' newline '\% complete' newline '$(\%c = 1 - p_{\mathrm{inc}}^{\mathrm{test}})$'], 'fontsize', 30); 
hide_only_axis('y');

subplot(3,1,[2,3]); hold on; 

% normalize for easy viewing, only for demo 
lgdn_objs(1) = plot(demo_res.t, normalize_minmax(demo_res.a), '-k', 'Displayname', '$a$');
lgdn_objs(2) = plot(demo_res.t, normalize_minmax(demo_res.bar_alpha_W), ':k', 'Displayname', '$\overline{\alpha_W}$'); 

xline(train_transition, ':', 'linewidth', 3);
scatter(select_steps, -0.5*ones(length(select_steps),1), 500, cmap, 'filled');

ylim([-0.5,1.1]);

yyaxis right; 
% smooth for demo 
lgdn_objs(3) = plot(demo_res.t, demo_res.b, '-', 'Color', 0.7*[1,1,1], 'Displayname', '$b$'); 
set(gca, 'visible', 'off');

lgdn_objs = lgdn_objs([1,3,2]);
lgnd_handle = legend(lgdn_objs, 'interpreter', 'latex', 'fontsize', 40, 'location', 'westoutside', 'numcolumns', 1);
lgnd_handle.Position = lgnd_handle.Position + [-0.08,0,0,0];

linkaxes(findall(gcf,'type','axes'),'x');
xlim([2e3,3.1e5]);

fig_name = fullfile(fig_path, fig_prefix);
export_fig(fig_name, '-r300', '-p0.02');
close;

%% Demo activation function 
fig_prefix = 'demo-training-select-actfun';
x = -0.5:0.01:2;
y = arrayfun(@(a,b) sigmoidal_activation(x,a,b), a_select, b_select, 'uni', 0); 
y = vertcat(y{:});

figure('WindowStyle', 'normal', 'DefaultLineLineWidth', 6, ...
    'Units', 'normalized', 'Position', [0.1,0.1,0.5,0.5]); 
hold on;
colororder(gca, cmap); 
plot(x,y); 
set(gca, 'xcolor', 'none', 'ycolor', 'none');
title('activation function', 'fontsize', 50);

fig_name = fullfile(fig_path, fig_prefix);
export_fig(fig_name, '-r300', '-p0.02');
close;