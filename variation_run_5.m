clc; clear; close all; 
run start_up; 

script_filename = mfilename;
file_prefix = 'full_var-sigmoidwiththres-IP-onoff_var-init_var-psigmoid_IPfasterthanSP'; 
max_eta_ip_a = 1e-2;
max_eta_ip_theta = 1e-3; 

%% Default opts 
def_opts = struct; 
def_opts.N = 10;

def_opts.mu = 1/(2*def_opts.N);
def_opts.selrow = 8; 

def_opts.train_norm_input = true;  % true or false
% def_opts.num_train = 30e4; 
% def_opts.p_train_complete = 1/15;
def_opts.num_train = 20e4; 
def_opts.p_train_complete = 1/10;

def_opts.eta_sp_ltp = 1e-3; 
def_opts.kappa_sp = 1/20; % kappa_sp=eta_sp_ltd/eta_sp_ltp = 0, 1/20, 1/15, 1/10

def_opts.weight_properties.norm_order = 2; % 1 or 2
def_opts.weight_properties.onlyexcitatory = true; % true or false 

def_opts.subsampled = 1000;

def_opts.num_test_per_pinc = 1000;
def_opts.test_p_sel = 0.5; % 50-50 of selective input 
def_opts.test_sigma_noise = 0.1; % small gaussian noise 
% def_opts.test_p_incs = [0, 0.1, 0.3, 0.5]; 
def_opts.test_p_incs = [0, 0.2, 0.4]; 

def_opts.test_shuffle = false; % not necessary 

def_opts.train_norm_input = true;
def_opts.test_norm_input  = true;

%% Variations 
num_rand = 10;  
% train_p_inc_vec = [0,0.1,0.3,0.5];
train_p_inc_vec = [0,0.2,0.4];

eta_ip_a_vec = [0, max_eta_ip_a];
eta_ip_theta_vec = [0,max_eta_ip_theta];

a_init_vec = linspace(1,10,10);
theta_init_vec = linspace(0.1,0.5,10);

p_sigmoid_vec = [0.2]; 

params_to_vary = {...
    struct('label', 'eta_ip_a', 'vec', eta_ip_a_vec), ...
    struct('label', 'eta_ip_theta', 'vec', eta_ip_theta_vec), ...
    struct('label', 'a_init', 'vec', a_init_vec), ...
    struct('label', 'theta_init', 'vec', theta_init_vec), ...
    struct('label', 'p_sigmoid', 'vec', p_sigmoid_vec), ...
    struct('label', 'train_p_inc', 'vec', train_p_inc_vec), ...
    };

[param_combs, param_labels, comb_indices] = return_combination(params_to_vary{:});

num_combs = size(param_combs,1);
num_params = length(param_labels);


%% Variables to save 
vars_tosave = {'a','theta', 'bar_alpha_W', 'dY_objvsnoise_test'};

%% Train 
train_results = cell(num_combs, 1); 
% ppm = ParforProgMon('Monitor', num_combs); 
ppm = ParforProgressbar(num_combs, ...
    'progressBarUpdatePeriod', 20);
    
t0_glob = tic;
parfor i = 1:num_combs
    t0i = tic; 
    opts = def_opts;
    for k = 1:num_params
        opts.(param_labels{k}) = param_combs(i,k);  %#ok<PFBNS>
    end
    
    result_array = arrayfun(@(~) train_with_incomplete_sigmoidthres(opts), 1:num_rand, 'uni', 1);
    result_aggregate = struct;
    
    for k = 1:length(vars_tosave)
        cat_vec = cat(1,result_array.(vars_tosave{k}));
        result_aggregate.(vars_tosave{k}).mean = squeeze(mean(cat_vec));
        result_aggregate.(vars_tosave{k}).sem = squeeze(std(cat_vec) / sqrt(num_rand));
    end
    
    result_aggregate.opts = opts;
    train_results{i} = result_aggregate;
    
    ppm.increment(); %#ok<PFBNS>
end

toc(t0_glob);

train_results = vertcat(train_results{:}); 

delete(ppm);

%% Save data 
save(fullfile('data', file_prefix), ...
    'script_filename', 'def_opts', 'max_eta_ip*', 'train_results', ...
    'param_combs', 'param_labels', 'comb_indices', 'params_to_vary', ...
    'vars_tosave');
