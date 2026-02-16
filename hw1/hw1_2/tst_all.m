
% Compare all three methods: Original CN, Permuted CN, ADI
% for accuracy and timing

fprintf('\n===== Original CN (Cholesky, no permutation) =====\n');
e2=1;
t0=tic;
N=80;  dt=.16000; heat2d
N=80;  dt=.08000; heat2d
N=80;  dt=.04000; heat2d
N=80;  dt=.02000; heat2d
N=80;  dt=.01000; heat2d
N=80;  dt=.00500; heat2d
N=80;  dt=.00250; heat2d
N=80;  dt=.00125; heat2d
etime_orig=toc(t0);
fprintf('Total time (original): %f s\n\n', etime_orig);

fprintf('===== Permuted CN (symamd + Cholesky) =====\n');
e2=1;
t0=tic;
N=80;  dt=.16000; heat2d_perm
N=80;  dt=.08000; heat2d_perm
N=80;  dt=.04000; heat2d_perm
N=80;  dt=.02000; heat2d_perm
N=80;  dt=.01000; heat2d_perm
N=80;  dt=.00500; heat2d_perm
N=80;  dt=.00250; heat2d_perm
N=80;  dt=.00125; heat2d_perm
etime_perm=toc(t0);
fprintf('Total time (permuted): %f s\n\n', etime_perm);

fprintf('===== ADI =====\n');
e2=1;
t0=tic;
N=80;  dt=.16000; heat2d_adi
N=80;  dt=.08000; heat2d_adi
N=80;  dt=.04000; heat2d_adi
N=80;  dt=.02000; heat2d_adi
N=80;  dt=.01000; heat2d_adi
N=80;  dt=.00500; heat2d_adi
N=80;  dt=.00250; heat2d_adi
N=80;  dt=.00125; heat2d_adi
etime_adi=toc(t0);
fprintf('Total time (ADI):      %f s\n\n', etime_adi);

fprintf('===== Timing Summary =====\n');
fprintf('Original:  %f s\n', etime_orig);
fprintf('Permuted:  %f s  (speedup: %.2fx)\n', etime_perm, etime_orig/etime_perm);
fprintf('ADI:       %f s  (speedup: %.2fx)\n', etime_adi, etime_orig/etime_adi);
