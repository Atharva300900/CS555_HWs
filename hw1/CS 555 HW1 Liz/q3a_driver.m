% q3a_driver.m
% Q3a driver: compare baseline CN vs symamd-permuted CN
% Fixed dt, vary N. Reports errors, runtimes, and speedup.

clear all
clc

hdr

dt = 0.01;     % keep dt fixed
T  = 1.2;      % match reference heat2d.m final time
Ns = [40 80 120 160 200];

fprintf('\n===== Q3a: CN vs Permuted CN (dt=%.4f, T=%.3f) =====\n', dt, T);
fprintf('   N   n=(N-1)^2 |  err_CN    err_perm |  time_CN   time_perm | speedup\n');
fprintf('  ---  --------- | ---------  -------- | ---------  ---------  | -------\n');

for N = Ns
    n = (N-1)^2;

    t0 = tic;
    e_cn = heat2d_cn(N, dt, T);
    t_cn = toc(t0);

    t0 = tic;
    e_pm = heat2d_cn_perm(N, dt, T);
    t_pm = toc(t0);

    sp = t_cn / t_pm;

    fprintf(' %3d  %9d | %9.2e  %9.2e | %9.4f  %9.4f | %7.2fx\n', ...
            N, n, e_cn, e_pm, t_cn, t_pm, sp);
end