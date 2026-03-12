function cfg = hw2_config()
% Shared constants for CS555 HW2.
cfg.nu = 1 / (100 * pi);
cfg.T = 2.0;
cfg.reference_s_star = 152.00516;
cfg.reference_t_star_paper = 1.6037;
cfg.reference_t_star = cfg.reference_t_star_paper / pi;
cfg.fixed_cfl = 0.1;
cfg.N_q1a = 200;
cfg.N_list = 2 .^ (5:11);
end
