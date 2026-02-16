% Plot for 2c: Convergence comparison ADI-CN vs ADI-BDF3

hdr

dts  = [0.16 0.08 0.04 0.02 0.01 0.005 0.0025 0.00125];
e_cn   = zeros(size(dts));
e_bdf3 = zeros(size(dts));

for idx = 1:length(dts)
  dt = dts(idx);
  clear N; N = 80;
  e2=1; heat2d_adi;      e_cn(idx) = e2;
  clear N; N = 80;
  e2=1; heat2d_adi_bdf3; e_bdf3(idx) = e2;
end

figure(1);
loglog(dts, e_cn,   'ro-', 'linewidth', 2, 'markersize', 8); hold on;
loglog(dts, e_bdf3, 'bs-', 'linewidth', 2, 'markersize', 8);

% Reference slopes
loglog(dts, 0.5*dts.^2, 'r--', 'linewidth', 1);
loglog(dts, 2.0*dts.^3, 'b--', 'linewidth', 1);
hold off;

xlabel('\Delta t','FontSize',14);
ylabel('Relative L_2 error','FontSize',14);
title('Convergence: ADI-CN (2nd order) vs ADI-BDF3 (3rd order)','FontSize',14);
legend('ADI-CN','ADI-BDF3','O(\Delta t^2) ref','O(\Delta t^3) ref','Location','SouthEast');
legend boxoff;
grid on;

print('-dpng','-r150','plot_q2c.png');
