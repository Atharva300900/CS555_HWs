
% Test ADI-BDF3 convergence and compare to ADI-CN

fprintf('\n===== ADI-CN (2nd order) =====\n');
e2=1; N=80;
dt=.16000; heat2d_adi
dt=.08000; heat2d_adi
dt=.04000; heat2d_adi
dt=.02000; heat2d_adi
dt=.01000; heat2d_adi
dt=.00500; heat2d_adi
dt=.00250; heat2d_adi
dt=.00125; heat2d_adi

fprintf('\n===== ADI-BDF3 (3rd order) =====\n');
e2=1; N=80;
dt=.16000; heat2d_adi_bdf3
dt=.08000; heat2d_adi_bdf3
dt=.04000; heat2d_adi_bdf3
dt=.02000; heat2d_adi_bdf3
dt=.01000; heat2d_adi_bdf3
dt=.00500; heat2d_adi_bdf3
dt=.00250; heat2d_adi_bdf3
dt=.00125; heat2d_adi_bdf3

% Timing comparison
fprintf('\n===== Timing Comparison =====\n');
fprintf('  N       ADI-CN      ADI-BDF3\n');
for Nval = [40 80 120 160 200]
  dt = 0.01;
  clear N; N = Nval;
  e2=1; t0=tic; heat2d_adi;      t_cn=toc(t0);
  clear N; N = Nval;
  e2=1; t0=tic; heat2d_adi_bdf3; t_bdf3=toc(t0);
  fprintf('  %3d   %8.4f s   %8.4f s\n', Nval, t_cn, t_bdf3);
end
