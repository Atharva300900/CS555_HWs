
% Scaling test: call heat2d, heat2d_perm, and heat2d_adi for each N, compare timings

dt = 0.01;

fprintf('\n===== Timing vs N: Original CN vs Permuted CN vs ADI (dt=%.4f) =====\n', dt);
fprintf('  N    n=(N-1)^2  |  orig(s)   perm(s)   adi(s)   spdup_perm  spdup_adi\n');
fprintf('  ----  ---------- | ---------  --------  --------  ----------  ---------\n');

e2 = 1;
for N = [40,80,120,160,200]

  t0 = tic; heat2d;      t_orig = toc(t0);
  t0 = tic; heat2d_perm; t_perm = toc(t0);
  t0 = tic; heat2d_adi;  t_adi  = toc(t0);
  t0 = tic; [e2, t_fac_out] = heat2d_adi_2(N, dt, 1.2);  t_adi_2  = toc(t0);

  fprintf('  %3d   %6d     |  %7.3fs   %7.3fs   %7.3fs   %7.3fs %8.1fx   %8.1fx\n', ...
          N, (N-1)^2, t_orig, t_perm, t_adi, t_adi_2, t_orig/t_adi, t_orig/t_adi_2);
end
