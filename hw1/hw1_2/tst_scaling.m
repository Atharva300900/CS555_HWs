
% Scaling test: compare timings as N increases

fprintf('\n===== Scaling Comparison (dt=0.01) =====\n');
fprintf('  N       n=N^2   Original    Permuted      ADI       Sp(Perm) Sp(ADI)\n');

for Nval = [40 80 120 160 200]
  dt = 0.01;

  clear N; N = Nval;
  e2=1; t0=tic; heat2d;      t_orig=toc(t0);

  clear N; N = Nval;
  e2=1; t0=tic; heat2d_perm; t_perm=toc(t0);

  clear N; N = Nval;
  e2=1; t0=tic; heat2d_adi;  t_adi=toc(t0);

  fprintf('  %3d   %6d   %8.3f s   %8.3f s   %8.4f s    %5.1fx    %6.1fx\n', ...
          Nval, (Nval-1)^2, t_orig, t_perm, t_adi, t_orig/t_perm, t_orig/t_adi);
end
