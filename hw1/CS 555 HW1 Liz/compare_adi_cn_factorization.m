% compare_adi_cn_factorization.m
% Compare ADI-CN with Cholesky vs without explicit factorization

clear all
clc

hdr

dt = 0.01;
T  = 1.2;
Ns = [40 80 120 160 200];

fprintf('\n===== ADI-CN: with vs without Cholesky (dt=%.4f, T=%.3f) =====\n', dt, T);
fprintf('   N   n=(N-1)^2 | time_chol  time_nofac | speedup\n');
fprintf('  ---  --------- | ---------- ----------- | --------\n');

for N = Ns
    n = (N-1)^2;

    t0 = tic;
    heat2d_adi_cn(N, dt, T);
    t_chol = toc(t0);

    t0 = tic;
    heat2d_adi_cn_nofac(N, dt, T);
    t_nofac = toc(t0);

    sp = t_nofac / t_chol;

    fprintf(' %3d  %9d | %10.4f %11.4f | %7.2fx\n', ...
            N, n, t_chol, t_nofac, sp);
end