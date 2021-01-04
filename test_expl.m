clc; clear; %close all; 

run start_up.m

fig_path = 'figures/misc/explanation/';

%%
mu_Y = 0.1;
eta_ip_a = 1e-3;
eta_ip_b = 1e-3; 

T = 2e4; 

N = 5;

t_switch = T/2; 
% inp_fun = @(N,t) 5+(t>t_switch)*linspace(0,-4,N)'; 
inp_fun = @(N,t) 5+randn(N,1)+(t>t_switch)*linspace(0,-4,N)'; 
% inp_fun = @(N,t) ((t>t_switch)*linspace(0,5,N)' + 1).*randn(N,1);

% fig_name = 'single-cell-noip-nonoiseinp';
% fig_name = 'single-cell-gainplast-nonoiseinp';
% fig_name = 'single-cell-biasplast-nonoiseinp';
% fig_name = 'single-cell-fullip-nonoiseinp';

% fig_name = 'single-cell-noip-noiseinp';
% fig_name = 'single-cell-gainplast-noiseinp';
% fig_name = 'single-cell-biasplast-noiseinp';
% fig_name = 'single-cell-fullip-noiseinp';


state_monitor.x = zeros(N,T);
state_monitor.y = zeros(N,T); 
state_monitor.a = zeros(N,T);
state_monitor.b = zeros(N,T); 
state_monitor.c = zeros(N,T); 

a = 1;
b = -5; 
for t = 1:T
    x = inp_fun(N,t); 
    y = sigmoidal_activation(x, a, b);
    
    common_factor = 1 - (2+1./mu_Y).*y + (y.^2)./mu_Y;
    da = eta_ip_a * (1./a + x .* common_factor);
    db = eta_ip_b * common_factor;
    a = a + da;
    b = b + db;
    
    
    state_monitor.x(:,t) = x;
    state_monitor.y(:,t) = y;
    state_monitor.a(:,t) = a;
    state_monitor.b(:,t) = b;
    state_monitor.c(:,t) = (common_factor);
    
    
end


smooth_opts.x = 10;
smooth_opts.y = 10;
smooth_opts.a = 0;
smooth_opts.b = 0;
smooth_opts.c = 100;

t = 1:T; 
cmap = parula(N)*0.9; 
state_vars = fieldnames(state_monitor); 
num_statevars = length(state_vars);
nrows = 3;
ncols = ceil(num_statevars/nrows); 

figure; 
for i = 1:num_statevars
    var_name = state_vars{i}; 
    var_vec  = state_monitor.(var_name); 
    if smooth_opts.(var_name) > 0
        var_vec = smoothdata(var_vec, 2, 'gaussian', smooth_opts.(var_name));
    end
    
    subplot(nrows, ncols, i); hold on; 
    colororder(gca, cmap); 
    plot(t, var_vec);
    title(sprintf('$%s$', var_name)); 
    colormap(cmap);
    colorbar; 
end

linkaxes(findall(gcf, 'type', 'axes'), 'x');
xlim([0.45,0.6]*T);

subplot(nrows,ncols,nrows*ncols); hold on; 
colororder(gca, cmap);
arrayfun(@(i) histogram(state_monitor.y(i,ceil(T*4/5):T), 'Normalization', 'pdf'), 1:N);
set(gca, 'yscale', 'log');

% export_fig(fullfile(fig_path, fig_name), '-dpng', '-r200', '-p0.02')
% close;

%%
% a = 0.5; 
% b = -4; 
% 
% x = 0:0.01:10;
% y_1 = sigmoidal_activation(x,a,b); 
% y_2 = sigmoidal_activation(x,a+1,b); 
% c_1 = x .* ( 1 - (2+1./mu_Y).*y_1 + (y_1.^2)./mu_Y );
% c_2 = x .* ( 1 - (2+1./mu_Y).*y_2 + (y_2.^2)./mu_Y );
% 
% figure; hold on; 
% plot(x, y_1, '-b');
% plot(x, y_2, '-r');

%%
a = 6; 
b = -4.5; 

x = -0.5:0.01:2; 
y_0 = sigmoidal_activation(x, a, b);
y_a = sigmoidal_activation(x, a+1, b);
y_b = sigmoidal_activation(x, a, b+1);

figure('units', 'normalized', 'position', [0,0,0.5,0.8])
subplot(211); hold on;
plot(x, y_0, '-k');
plot(x, y_a, '-r');
plot(x, y_b, '-b');
title('act fun');
legend('base', 'inc a', 'inc b');

subplot(212); hold on;
plot(x, y_0 - y_0, '-k');
plot(x, y_a - y_0, '-r');
plot(x, y_b - y_0, '-b');
title('subtracted by base');
legend('base', 'inc a', 'inc b');

fig_name = 'actfun-subbase';
export_fig(fullfile(fig_path, fig_name), '-dpng', '-r200', '-p0.02')
close;

%%
x = -0.5:0.01:2; 
y_0 = sigmoidal_activation(x, a, b);
y_a = sigmoidal_activation(x, a+5, b);
y_b = sigmoidal_activation(x, a, b+1);

figure('units', 'normalized', 'position', [0,0,0.4,0.8])
subplot(511); hold on;
plot(x, y_0, '-k');
plot(x, y_a, '-r');
plot(x, y_b, '-b');
title('act fun'); set(gca, 'xcolor', 'none')
legend('base', 'inc a by 5', 'inc b by 1');

subplot(5,1,[2,3]); hold on;
plot(x, y_0 - y_0, '-k');
plot(x, y_a - y_0, '-r');
plot(x, y_b - y_0, '-b');
title('subtracted by base'); set(gca, 'xcolor', 'none')
legend('base', 'inc a by 5', 'inc b by 1');

subplot(5,1,[4,5]); hold on;
plot(x, y_0 ./ y_0, '-k');
plot(x, y_a ./ y_0, '-r');
plot(x, y_b ./ y_0, '-b');
title('divided by base'); 
legend('base', 'inc a by 5', 'inc b by 1');

despline('all');

fig_name = 'actfun-normbase-inca5incb1';
export_fig(fullfile(fig_path, fig_name), '-dpng', '-r200', '-p0.02')
close;

%%
a = 6; 
b = -4.5; 
p_x = 0.7;
eta = 1e-2; 
T = 1e4; 

res = test_2x1y(p_x, a, b, eta, T);

figure; 
subplot(211); hold on;
plot(res.alpha_L1, '-');
plot(res.alpha_L2, ':');

subplot(212); hold on;
plot(res.w1, '-r');
plot(res.w2, '-k');

%%
a = 6; 
b = -4.5; 
p_x = 0.6;
eta = 1e-3; 
T = 3e4; 
n = 10;

var_s = struct;

var_s.base.a = a;
var_s.base.b = b;
var_s.base.color = 'k';
var_s.base.name = 'base';

var_s.GAIN.a = a+1;
var_s.GAIN.b = b;
var_s.GAIN.color = 'r';
var_s.GAIN.name = 'a \leftarrow a+1';

var_s.BIAS.a = a;
var_s.BIAS.b = b+1;
var_s.BIAS.color = 'b';
var_s.BIAS.name = 'b \leftarrow b+1';

res_s = structfun(@(s) ...
        structfun(@(x) ...
            struct('mean', mean(cell2mat(x),2), 'sem', 2.58*std(cell2mat(x),0,2)/sqrt(n)), ...
            structarray_to_struct(...
                arrayfun(@(~) ...
                    test_2x1y(p_x, s.a, s.b, eta, T), ...
                1:n),...
            0), ...
        'uni', 0), ...
    var_s, 'uni', 0);

conds = fieldnames(res_s); 
for i = 1:length(conds)
    res_s.(conds{i}) = mergestruct(res_s.(conds{i}), var_s.(conds{i}));
end

%%

t = 1:T;

figure('units', 'normalized', 'position', [0,0,0.7,1])

subplot(211); hold on;
lgnd_objs_1 = structfun(@(r) plot(t, r.alpha_L1.mean, r.color, 'displayname', ...
    ['$\alpha_{L1} - ' r.name '$']), res_s);
structfun(@(r) fill([t fliplr(t)], [r.alpha_L1.mean-r.alpha_L1.sem;flipud(r.alpha_L1.mean+r.alpha_L1.sem)]', ...
    r.color, 'FaceAlpha', 0.2, 'LineStyle', 'none'), res_s);

lgnd_objs_2 = structfun(@(r) plot(t, r.alpha_L2.mean, r.color, 'linestyle', ':', ...
    'displayname',  ['$\alpha_{L2} - ' r.name '$']), res_s);
structfun(@(r) fill([t fliplr(t)], [r.alpha_L2.mean-r.alpha_L2.sem;flipud(r.alpha_L2.mean+r.alpha_L2.sem)]', ...
    r.color, 'FaceAlpha', 0.2, 'LineStyle', 'none'), res_s);
legend([lgnd_objs_1;lgnd_objs_2], 'Interpreter', 'latex', 'NumColumns', 2);


subplot(212); hold on;
lgnd_objs_1 = structfun(@(r) plot(t,r.w1.mean, '-', 'color', r.color, ...
     'displayname',  ['$w_1 - ' r.name '$']), res_s);
lgnd_objs_2 = structfun(@(r) plot(t,r.w2.mean, '--', 'color', r.color, ...
     'displayname',  ['$w_2 - ' r.name '$']), res_s);
structfun(@(r) fill([t fliplr(t)], [r.w1.mean-r.w1.sem;flipud(r.w1.mean+r.w1.sem)]', ...
    r.color, 'FaceAlpha', 0.2, 'LineStyle', 'none'), res_s);
structfun(@(r) fill([t fliplr(t)], [r.w2.mean-r.w2.sem;flipud(r.w2.mean+r.w2.sem)]', ...
    r.color, 'FaceAlpha', 0.2, 'LineStyle', 'none'), res_s);
legend([lgnd_objs_1;lgnd_objs_2], 'Interpreter', 'latex', ...
    'NumColumns', 2, 'Location', 'east');

despline('all');

fig_name = 'demo-2inps';
export_fig(fullfile(fig_path, fig_name), '-dpng', '-r200', '-p0.02')
close;

%%
x = 0:0.01:1;
alpha = 0.56; 
eta = 1; 

a = 6; 
b = -4.5; 
% a = 1; 
% b = 0; 

r0 = testfun(alpha,a,b,eta,x);
ra = testfun(alpha,a+1,b,eta,x);
rb = testfun(alpha,a,b+1,eta,x);

figure; hold on;
plot(x,r0,'-k')
plot(x,ra,'-r')
plot(x,rb,'-b')
xline(alpha);
yline(1);