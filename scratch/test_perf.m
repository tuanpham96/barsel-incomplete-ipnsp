opts = struct; 
opts.N = 10;
opts.a_init = 6; % 1 when norm_input = false, 5.7 otw
opts.b_init = -4.7; % -4.6 when norm_input = false, -4.7 otw
opts.mu = 1/(2*opts.N);
opts.selrow = 8; 

opts.train_norm_input = true;  % true or false
opts.num_train = 20e4; 
opts.p_train_complete = 1/2;

opts.num_test_per_pinc = 1000;
opts.test_p_sel = 0.5;
opts.test_sigma_noise = 0.1; 
opts.test_norm_input = true; 
opts.test_p_incs = [0, 0.2, 0.4, 0.6]; 
opts.test_shuffle = false; % not necessary 

opts.eta_sp_ltp = 0e-3; 
opts.kappa_sp = 0; % kappa_sp=eta_sp_ltd/eta_sp_ltp = 0, 1/20, 1/15, 1/10

opts.weight_properties.norm_order = 2; % 1 or 2
opts.weight_properties.onlyexcitatory = true; % true or false 

opts.eta_ip_b = 0e-3;
opts.subsampled = 1000;

opts.eta_ip_a = 1e-3;
opts.k_r = 0.4;
opts.train_p_inc = 0.1;
tic
res = train_with_incomplete(opts); 
toc

figure;

subplot(311); hold on; 
plot(res.t, res.a);

subplot(312); hold on; 
plot(res.t, res.b);

subplot(313); hold on; 
colororder(parula(4));
plot(res.t, squeeze(res.dY_objvsnoise_test));
colorbar;
%%
colororder(parula(4));
plot(res.t, squeeze(res.dY_objvsnoise_test), '--');
