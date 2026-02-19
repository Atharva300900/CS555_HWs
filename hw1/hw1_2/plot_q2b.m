% Plot for 2b: Timing scaling comparison (all three methods)

hdr

Nvals = [40 60 80 100 120 140 160 180 200];
t_orig = zeros(size(Nvals));
t_perm = zeros(size(Nvals));
t_adi  = zeros(size(Nvals));

for idx = 1:length(Nvals)
  dt = 0.01;
  clear N; N = Nvals(idx);
  e2=1; t0=tic; heat2d;      t_orig(idx)=toc(t0);
  clear N; N = Nvals(idx);
  e2=1; t0=tic; heat2d_perm; t_perm(idx)=toc(t0);
  clear N; N = Nvals(idx);
  e2=1; t0=tic; [e2, t_fac_out] = heat2d_adi_2(N, dt, 1.2);  t_adi(idx)=toc(t0);
end

figure(1);
semilogy(Nvals, t_orig, 'ko-', 'linewidth', 2, 'markersize', 8); hold on;
semilogy(Nvals, t_perm, 'rs-', 'linewidth', 2, 'markersize', 8);
semilogy(Nvals, t_adi,  'b^-', 'linewidth', 2, 'markersize', 8);
hold off;
xlabel('N','FontSize',14);
ylabel('Wall-clock time (s)','FontSize',14);
title('Timing Scaling: Original vs Permuted vs ADI','FontSize',14);
legend('Original CN','Permuted CN','ADI','Location','NorthWest');
legend boxoff;
grid on;

print('-dpng','-r150','plot_q2b_timing.png');
