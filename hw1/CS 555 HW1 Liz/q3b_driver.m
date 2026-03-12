% q3b_driver.m
% Q3b driver: compare CN, permuted CN, and ADI-CN
% Fixed dt, vary N. Reports errors, runtimes, speedups, and scalability plot.

close all
clear all
clc

hdr

dt = 0.01;
T  = 1.2;
Ns = [40 80 120 160 200];

t_cn   = zeros(size(Ns));
t_perm = zeros(size(Ns));
t_adi  = zeros(size(Ns));

e_cn   = zeros(size(Ns));
e_perm = zeros(size(Ns));
e_adi  = zeros(size(Ns));

fprintf('\n===== Q3b: ADI-CN vs CN Baselines (dt=%.4f, T=%.3f) =====\n', dt, T);
fprintf('   N   n=(N-1)^2 |  err_CN    err_perm   err_ADI |  time_CN   time_perm  time_ADI | sp(perm) sp(ADI)\n');
fprintf('  ---  --------- | ---------  ---------  ------- | ---------  ---------  -------- | -------- -------\n');

for i = 1:length(Ns)
    N = Ns(i);
    n = (N-1)^2;

    t0 = tic;
    e_cn(i) = heat2d_cn(N, dt, T);
    t_cn(i) = toc(t0);

    t0 = tic;
    e_perm(i) = heat2d_cn_perm(N, dt, T);
    t_perm(i) = toc(t0);

    t0 = tic;
    e_adi(i) = heat2d_adi_cn(N, dt, T);
    t_adi(i) = toc(t0);

    sp_perm = t_cn(i) / t_perm(i);
    sp_adi  = t_cn(i) / t_adi(i);

    fprintf(' %3d  %9d | %9.2e  %9.2e  %9.2e | %9.4f  %9.4f  %8.4f | %8.2fx %7.2fx\n', ...
            N, n, e_cn(i), e_perm(i), e_adi(i), ...
            t_cn(i), t_perm(i), t_adi(i), sp_perm, sp_adi);
end

% -----------------------------
% Scalability plot
% -----------------------------
figure('Position', [100 100 600 600]);
loglog(Ns, t_cn,   'o-', 'Color','k', 'MarkerSize',10, lw, 3); hold on
loglog(Ns, t_perm, 's-', 'Color','r', 'MarkerSize',10, lw, 3);
loglog(Ns, t_adi,  '^-', 'Color','b', 'MarkerSize',10, lw, 3);
grid on

xlabel('N', fs, 20)
ylabel('Runtime (seconds)', fs, 20)
legend('CN', 'CN + symamd', 'ADI-CN', 'Location', 'NorthWest', 'FontSize', 20)
title('Q3b: Runtime Scaling with N', fs, 30)
