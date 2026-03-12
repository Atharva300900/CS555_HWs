function results = solve_q1a()
% HW2 Q1a: plot u(x,t) for N = 200 at t = 0.0, 0.1, ..., 2.0.
setup_graphics();

nu = 1 / (100 * pi);
T = 2.0;
N = 200;
CFL = 0.1;
t_save = (0.0:0.1:2.0).';
root = fileparts(mfilename('fullpath'));

results.uniform = burgers_solve(N, nu, T, 'uniform', 'convective', CFL, t_save);
results.chebyshev = burgers_solve(N, nu, T, 'chebyshev', 'convective', CFL, t_save);

fig = figure('visible', 'off', 'Position', [100 100 1000 650]);
hold on;
cmap = jet(numel(t_save));
style_uniform = plot(nan, nan, 'k-', 'LineWidth', 1.4);
style_chebyshev = plot(nan, nan, 'k--', 'LineWidth', 1.4);

for k = 1:numel(t_save)
    color_k = cmap(k, :);
    plot(results.uniform.xb, results.uniform.u_snap(:, k), '-', 'Color', color_k, 'LineWidth', 1.0);
    plot(results.chebyshev.xb, results.chebyshev.u_snap(:, k), '--', 'Color', color_k, 'LineWidth', 1.0);
end

hold off;
grid on;
xlabel('x');
ylabel('u(x,t)');
title(sprintf('Q1a: N = %d, solid = uniform, dashed = Chebyshev', N));
xlim([0.0, 1.0]);
legend([style_uniform, style_chebyshev], {'Uniform grid', 'Chebyshev grid'}, 'Location', 'northwest');
colormap(cmap);
caxis([0.0, 2.0]);
cb = colorbar('eastoutside');
ylabel(cb, 'time t');
annotation('textbox', [0.13 0.01 0.55 0.05], 'String', 'All curves are plotted for t = 0.0, 0.1, ..., 2.0', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

results.figure_path = fullfile(root, 'q1a_solution.png');
print(fig, results.figure_path, '-dpng', '-r200');
close(fig);

fig_zoom = figure('visible', 'off', 'Position', [100 100 1000 650]);
hold on;
style_uniform_zoom = plot(nan, nan, 'k-', 'LineWidth', 1.4);
style_chebyshev_zoom = plot(nan, nan, 'k--', 'LineWidth', 1.4);

for k = 1:numel(t_save)
    color_k = cmap(k, :);
    plot(results.uniform.xb, results.uniform.u_snap(:, k), '-', 'Color', color_k, 'LineWidth', 1.0);
    plot(results.chebyshev.xb, results.chebyshev.u_snap(:, k), '--', 'Color', color_k, 'LineWidth', 1.0);
end

hold off;
grid on;
xlabel('x');
ylabel('u(x,t)');
title(sprintf('Q1a Zoom: N = %d, x in [0.990, 1.000]', N));
xlim([0.99, 1.0]);
legend([style_uniform_zoom, style_chebyshev_zoom], {'Uniform grid', 'Chebyshev grid'}, 'Location', 'southwest');
colormap(cmap);
caxis([0.0, 2.0]);
cb = colorbar('eastoutside');
ylabel(cb, 'time t');
annotation('textbox', [0.13 0.01 0.55 0.05], 'String', 'All curves are plotted for t = 0.0, 0.1, ..., 2.0', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

results.zoom_figure_path = fullfile(root, 'q1a_solution_zoom_x0990_to_1.png');
print(fig_zoom, results.zoom_figure_path, '-dpng', '-r200');
close(fig_zoom);

% Compare the two solutions on the uniform grid to visualize the boundary-layer discrepancy.
x_compare = results.uniform.xb;
u_chebyshev_on_uniform = zeros(size(results.uniform.u_snap));
for k = 1:numel(t_save)
    u_chebyshev_on_uniform(:, k) = interp1(results.chebyshev.xb, results.chebyshev.u_snap(:, k), ...
        x_compare, 'pchip');
end
u_diff = results.uniform.u_snap - u_chebyshev_on_uniform;

fig_diff = figure('visible', 'off', 'Position', [100 100 1000 750]);
subplot(2, 1, 1);
hold on;
for k = 1:numel(t_save)
    plot(x_compare, u_diff(:, k), 'Color', cmap(k, :), 'LineWidth', 1.0);
end
plot([0.0, 1.0], [0.0, 0.0], 'k:', 'LineWidth', 1.0);
hold off;
grid on;
xlabel('x');
ylabel('u_{uniform} - u_{Chebyshev interp}');
title(sprintf('Q1a Difference: uniform - Chebyshev on uniform grid, N = %d', N));
xlim([0.0, 1.0]);

subplot(2, 1, 2);
hold on;
for k = 1:numel(t_save)
    plot(x_compare, u_diff(:, k), 'Color', cmap(k, :), 'LineWidth', 1.0);
end
plot([0.9, 1.0], [0.0, 0.0], 'k:', 'LineWidth', 1.0);
hold off;
grid on;
xlabel('x');
ylabel('u_{uniform} - u_{Chebyshev interp}');
title('Zoom near x = 1');
xlim([0.9, 1.0]);

colormap(cmap);
caxis([0.0, 2.0]);
cb = colorbar('eastoutside');
ylabel(cb, 'time t');
annotation('textbox', [0.13 0.01 0.75 0.05], 'String', ...
    'Chebyshev values are interpolated onto the uniform grid using pchip before differencing.', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left');

results.x_compare = x_compare;
results.u_chebyshev_on_uniform = u_chebyshev_on_uniform;
results.u_diff = u_diff;
results.diff_figure_path = fullfile(root, 'q1a_difference_uniform_minus_chebyshev.png');
print(fig_diff, results.diff_figure_path, '-dpng', '-r200');
close(fig_diff);
end

function setup_graphics()
setenv('GNUTERM', 'pngcairo');
set(0, 'defaultfigurevisible', 'off');
end
