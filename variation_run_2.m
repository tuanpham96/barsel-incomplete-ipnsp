clc; clear; close all; 
run start_up; 

% file_prefix = 'full_var-IP-onoff_var-kr-var-abinit_var-norminp_IPnoSP'; 
file_prefix = 'full_var-IP-onoff_var-kr-var-abinit_var-norminp_IPnoSP'; 

max_eta_ip_a = 1e-2;
max_eta_ip_b = 1e-2;

%% Default opts 
def_opts = struct; 
def_opts.N = 10;

def_opts.mu = 1/(2*def_opts.N);
def_opts.selrow = 8; 

def_opts.train_norm_input = true;  % true or false
def_opts.num_train = 10e4; 
def_opts.p_train_complete = 0;

def_opts.eta_sp_ltp = 0;
def_opts.kappa_sp = 1/20; % kappa_sp=eta_sp_ltd/eta_sp_ltp = 0, 1/20, 1/15, 1/10

def_opts.weight_properties.norm_order = 2; % 1 or 2
def_opts.weight_properties.onlyexcitatory = true; % true or false 

def_opts.subsampled = 1000;

def_opts.num_test_per_pinc = 1000;
def_opts.test_p_sel = 0.5; % 50-50 of selective input 
def_opts.test_sigma_noise = 0.1; % small gaussian noise 
def_opts.test_p_incs = [0, 0.5]; 
def_opts.test_shuffle = false; % not necessary 

%% Variations of parameters
num_rand = 5; 

train_p_inc_vec = [0,0.5];

eta_ip_a_vec = [0,max_eta_ip_a];
eta_ip_b_vec = [0,max_eta_ip_b];

k_r_vec = [0.4,0.1]; 
a_init_vec = linspace(1,10,10);
b_init_vec = linspace(-10,0,10);

norm_input_opt_vec = [1]; 

params_to_vary = {...
    struct('label', 'eta_ip_a', 'vec', eta_ip_a_vec), ...
    struct('label', 'eta_ip_b', 'vec', eta_ip_b_vec), ...
    struct('label', 'k_r', 'vec', k_r_vec), ...
    struct('label', 'a_init', 'vec', a_init_vec), ...
    struct('label', 'b_init', 'vec', b_init_vec), ...
    struct('label', 'norm_input_opt', 'vec', norm_input_opt_vec), ...
    struct('label', 'train_p_inc', 'vec', train_p_inc_vec), ...
    };

[param_combs, param_labels, comb_indices] = return_combination(params_to_vary{:});

num_combs = size(param_combs,1);
num_params = length(param_labels);

%% Variables to save 
vars_tosave = {'a','b', 'dY_objvsnoise_test'};

%% Train 
train_results = cell(num_combs, 1); 
ppm = ParforProgMon('Monitor', num_combs); 

tic
parfor i = 1:num_combs
    
    opts = def_opts;
    for k = 1:num_params
        opts.(param_labels{k}) = param_combs(i,k);  %#ok<PFBNS>
    end
    
    result_array = arrayfun(@(~) train_with_incomplete(opts), 1:num_rand, 'uni', 1);
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

toc

train_results = vertcat(train_results{:}); 

%% Save data 
save(fullfile('data', file_prefix), ...
    'def_opts', 'max_eta_ip*', 'train_results', ...
    'param_combs', 'param_labels', 'comb_indices', 'params_to_vary', ...
    'vars_tosave');
