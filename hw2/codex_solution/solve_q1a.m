function results = solve_q1a()
% HW2 Q1a: plot u(x,t) for N = 200 at t = 0.0, 0.1, ..., 2.0.
hw2_setup_graphics();
cfg = hw2_config();
paths = hw2_paths();

N = cfg.N_q1a;
target_dt = cfg.fixed_cfl / N;
t_save = (0.0:0.1:2.0).';
results.uniform = hw2_load_or_solve(N, 'uniform', 'convective', target_dt, t_save);
results.chebyshev = hw2_load_or_solve(N, 'chebyshev', 'convective', target_dt, t_save);

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

results.figure_path = fullfile(paths.figures, 'q1a_solution.png');
print(fig, results.figure_path, '-dpng', '-r200');
close(fig);

save(fullfile(paths.data, 'q1a_results.mat'), 'results');
end
