%% Question 1d: Conservation form comparison
% Test conservation form at finest resolution and compare with convective form

clear; close all;

nu = 1/(100*pi);
T = 2.0;
CFL = 0.1;
s_analytical = 152.00516;

N = 1600;  % finest resolution from 1c

fprintf('Question 1d: Conservation form vs Convective form at N=%d\n\n', N);

% Convective form
fprintf('Running convective form...\n');
res_conv = burgers_solve(N, nu, T, 'uniform', 'convective', CFL);
fprintf('  Convective: s* = %.5f, rel.err = %.4e\n', ...
    res_conv.s_star, abs(res_conv.s_star - s_analytical)/s_analytical);

% Conservation form
fprintf('Running conservation form...\n');
res_cons = burgers_solve(N, nu, T, 'uniform', 'conservation', CFL);
fprintf('  Conservation: s* = %.5f, rel.err = %.4e\n', ...
    res_cons.s_star, abs(res_cons.s_star - s_analytical)/s_analytical);

fprintf('\nDifference between forms: %.5e\n', abs(res_conv.s_star - res_cons.s_star));

% Also test with Chebyshev spacing
fprintf('\nWith Chebyshev spacing:\n');
res_conv_ch = burgers_solve(N, nu, T, 'chebyshev', 'convective', CFL);
fprintf('  Convective: s* = %.5f, rel.err = %.4e\n', ...
    res_conv_ch.s_star, abs(res_conv_ch.s_star - s_analytical)/s_analytical);

res_cons_ch = burgers_solve(N, nu, T, 'chebyshev', 'conservation', CFL);
fprintf('  Conservation: s* = %.5f, rel.err = %.4e\n', ...
    res_cons_ch.s_star, abs(res_cons_ch.s_star - s_analytical)/s_analytical);

% Plot comparison of s(t) for both forms
figure('Position', [100 100 800 500]);
plot(res_conv.t_all, res_conv.s_all, 'b-', 'LineWidth', 1.5); hold on;
plot(res_cons.t_all, res_cons.s_all, 'r--', 'LineWidth', 1.5);
xlabel('t', 'FontSize', 14);
ylabel('s(t)', 'FontSize', 14);
title(sprintf('Convective vs Conservation Form, N=%d, Uniform', N), 'FontSize', 14);
legend('Convective', 'Conservation', 'FontSize', 12);
grid on;
saveas(gcf, 'q1d_form_comparison.png');
fprintf('\nFigure saved to q1d_form_comparison.png\n');
