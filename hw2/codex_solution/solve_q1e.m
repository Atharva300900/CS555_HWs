function studies = solve_q1e()
% HW2 Q1e: convergence tables for four grid/form combinations.
hw2_setup_graphics();
cfg = hw2_config();
paths = hw2_paths();

case_defs = {
    'uniform',   'convective',   'q1e_uniform_convective';
    'chebyshev', 'convective',   'q1e_chebyshev_convective';
    'uniform',   'conservation', 'q1e_uniform_conservation';
    'chebyshev', 'conservation', 'q1e_chebyshev_conservation';
};

studies = struct();
fig = figure('visible', 'off', 'Position', [100 100 980 620], 'Color', 'w');
hold on;
legend_entries = cell(size(case_defs, 1), 1);
s_star_matrix = zeros(numel(cfg.N_list), size(case_defs, 1));
plot_colors = [
    0.00 0.27 0.58;
    0.75 0.18 0.12;
    0.10 0.48 0.22;
    0.45 0.20 0.58
];
plot_markers = {'o', 's', 'd', '^'};
plot_styles = {'-', '-', '--', '--'};

for c = 1:size(case_defs, 1)
    grid_name = case_defs{c, 1};
    form_name = case_defs{c, 2};
    tag_name = case_defs{c, 3};
    rows = cell(numel(cfg.N_list), 5);
    s_vals = zeros(numel(cfg.N_list), 1);
    err_vals = zeros(numel(cfg.N_list), 1);
    ratio_vals = nan(numel(cfg.N_list), 1);
    nsteps_vals = zeros(numel(cfg.N_list), 1);

    for i = 1:numel(cfg.N_list)
        N = cfg.N_list(i);
        target_dt = cfg.fixed_cfl / N;
        result = hw2_load_or_solve(N, grid_name, form_name, target_dt, []);
        s_vals(i) = result.s_star;
        s_star_matrix(i, c) = result.s_star;
        err_vals(i) = abs(result.s_star - cfg.reference_s_star) / cfg.reference_s_star;
        nsteps_vals(i) = result.nsteps;
        if i > 1 && err_vals(i) > 0.0
            ratio_vals(i) = err_vals(i - 1) / err_vals(i);
            ratio_str = sprintf('%.3f', ratio_vals(i));
        else
            ratio_str = '--';
        end
        rows(i, :) = {sprintf('%d', N), sprintf('%d', result.nsteps), sprintf('%.8f', result.s_star), ...
            sprintf('%.3e', err_vals(i)), ratio_str};
    end

    key = sprintf('%s_%s', grid_name, form_name);
    studies.(key).N = cfg.N_list;
    studies.(key).nsteps = nsteps_vals;
    studies.(key).s_star = s_vals;
    studies.(key).rel_err = err_vals;
    studies.(key).ratio = ratio_vals;

    hw2_write_table_txt(fullfile(paths.tables, [tag_name '.txt']), {'N', 'nsteps', 's_star', 'rel_err', 'ratio'}, rows);
    figure(fig);
    loglog(cfg.N_list, err_vals, ...
        'Color', plot_colors(c, :), ...
        'LineStyle', plot_styles{c}, ...
        'LineWidth', 2.2, ...
        'Marker', plot_markers{c}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', plot_colors(c, :), ...
        'MarkerEdgeColor', [0.1 0.1 0.1]);
    legend_entries{c} = sprintf('%s / %s', capitalize_word(grid_name), capitalize_word(form_name));
end

figure(fig);
hold off;
grid on;
grid minor;
xlim([min(cfg.N_list), max(cfg.N_list)]);
set(gca, ...
    'XColor', [0.1 0.1 0.1], ...
    'YColor', [0.1 0.1 0.1], ...
    'GridColor', [0.75 0.75 0.75], ...
    'MinorGridColor', [0.88 0.88 0.88], ...
    'GridAlpha', 0.5, ...
    'MinorGridAlpha', 0.5, ...
    'LineWidth', 1.0, ...
    'FontSize', 12, ...
    'Box', 'on', ...
    'Layer', 'top');
xlabel('N');
ylabel('relative error in s^*');
title('Q1e: convergence for the four required cases');
legend(legend_entries, 'Location', 'southwest', 'Box', 'off', 'FontSize', 11);
studies.figure_path = fullfile(paths.figures, 'q1e_convergence.png');
print(fig, studies.figure_path, '-dpng', '-r260');
close(fig);

fig_ref = figure('visible', 'off', 'Position', [100 100 980 620], 'Color', 'w');
hold on;
for c = 1:size(case_defs, 1)
    plot(cfg.N_list, s_star_matrix(:, c), ...
        'Color', plot_colors(c, :), ...
        'LineStyle', plot_styles{c}, ...
        'LineWidth', 2.2, ...
        'Marker', plot_markers{c}, ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', plot_colors(c, :), ...
        'MarkerEdgeColor', [0.1 0.1 0.1]);
end
plot(cfg.N_list, cfg.reference_s_star * ones(size(cfg.N_list)), 'k-.', 'LineWidth', 1.8);
hold off;
set(gca, 'XScale', 'log');
grid on;
grid minor;
set(gca, ...
    'XColor', [0.1 0.1 0.1], ...
    'YColor', [0.1 0.1 0.1], ...
    'GridColor', [0.75 0.75 0.75], ...
    'MinorGridColor', [0.88 0.88 0.88], ...
    'GridAlpha', 0.5, ...
    'MinorGridAlpha', 0.5, ...
    'LineWidth', 1.0, ...
    'FontSize', 12, ...
    'Box', 'on', ...
    'Layer', 'top');
xlabel('N');
ylabel('s^*');
title('Q1e: computed s^* versus analytical reference');
legend_entries_ref = legend_entries;
legend_entries_ref{end + 1} = 'Analytical s^*';
legend(legend_entries_ref, 'Location', 'southeast', 'Box', 'off', 'FontSize', 11);
studies.reference_figure_path = fullfile(paths.figures, 'q1e_s_star_vs_reference.png');
print(fig_ref, studies.reference_figure_path, '-dpng', '-r260');
close(fig_ref);

save(fullfile(paths.data, 'q1e_results.mat'), 'studies');
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end
