function summary = solve_q1d()
% HW2 Q1d: compare convective and conservation forms at the finest N.
hw2_setup_graphics();
cfg = hw2_config();
paths = hw2_paths();
N = cfg.N_list(end);
target_dt = cfg.fixed_cfl / N;
grids = {'uniform', 'chebyshev'};
summary = struct();
rows = cell(numel(grids), 5);

fig = figure('visible', 'off', 'Position', [100 100 950 650]);

for i = 1:numel(grids)
    grid_name = grids{i};
    convective = hw2_load_or_solve(N, grid_name, 'convective', target_dt, []);
    conservative = hw2_load_or_solve(N, grid_name, 'conservation', target_dt, []);

    summary.(grid_name).convective = convective;
    summary.(grid_name).conservative = conservative;
    summary.(grid_name).difference = abs(convective.s_star - conservative.s_star);

    rows(i, :) = {capitalize_word(grid_name), sprintf('%.8f', convective.s_star), sprintf('%.8f', conservative.s_star), ...
        sprintf('%.3e', abs(conservative.s_star - cfg.reference_s_star) / cfg.reference_s_star), ...
        sprintf('%.3e', summary.(grid_name).difference)};

    subplot(2, 1, i);
    plot(convective.t_all, convective.s_all, 'b-', 'LineWidth', 1.3); hold on;
    plot(conservative.t_all, conservative.s_all, 'r--', 'LineWidth', 1.3);
    plot([0.0, cfg.T], [cfg.reference_s_star, cfg.reference_s_star], 'k:', 'LineWidth', 1.0);
    hold off;
    grid on;
    xlabel('t');
    ylabel('s(t)');
    title(sprintf('Q1d: %s grid, N = %d', capitalize_word(grid_name), N));
    legend({'Convective', 'Conservation', 'Reference s^*'}, 'Location', 'southeast');
end

summary.figure_path = fullfile(paths.figures, 'q1d_form_comparison.png');
print(fig, summary.figure_path, '-dpng', '-r200');
close(fig);

summary.table_path = fullfile(paths.tables, 'q1d_summary.txt');
hw2_write_table_txt(summary.table_path, {'Grid', 's_star_convective', 's_star_conservation', 'rel_err_conservation', 'abs_difference'}, rows);
save(fullfile(paths.data, 'q1d_summary.mat'), 'summary');
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end
