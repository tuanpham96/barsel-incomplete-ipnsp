opts = struct; 
opts.N = 10;
opts.a_init = 5.7; % 1 when norm_input = false, 5.7 otw
opts.b_init = -4.7; % -4.6 when norm_input = false, -4.7 otw
opts.mu = 1/(2*opts.N);
opts.selrow = 8; 

opts.train_norm_input = false;  % true or false
opts.num_train = 20e4; 
opts.p_train_complete = 1/10;

opts.num_test_per_pinc = 1000;
opts.test_p_sel = 0.5;
opts.test_sigma_noise = 0.1; 
opts.test_norm_input = false; 
opts.test_p_incs = [0, 0.2, 0.4, 0.6]; 
opts.test_shuffle = false; % not necessary 

opts.eta_sp_ltp = 1e-3; 
opts.kappa_sp = 1/20; % kappa_sp=eta_sp_ltd/eta_sp_ltp = 0, 1/20, 1/15, 1/10

opts.weight_properties.norm_order = 2; % 1 or 2
opts.weight_properties.onlyexcitatory = true; % true or false 

opts.eta_ip_b = 0e-2;
opts.subsampled = 1000;

opts.p_inc = 0.5;
opts.eta_ip_a = opts.eta_ip_b;

tic
res = train_with_incomplete(opts); 
toc
%%
figure; hold on; 
colororder(parula(4));
plot(res.t, squeeze(res.dY_objvsnoise_test));
colorbar;
%%
colororder(parula(4));
plot(res.t, squeeze(res.dY_objvsnoise_test), '--');
