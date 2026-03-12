function studies = solve_q1e()
% HW2 Q1e: convergence tables for four grid/form combinations.
setup_graphics();

nu = 1 / (100 * pi);
T = 2.0;
CFL = 0.1;
reference_s_star = 152.00516;
N_list = 2 .^ (5:11);
root = fileparts(mfilename('fullpath'));

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
    rows = cell(numel(N_list), 5);
    s_vals = zeros(numel(N_list), 1);
    err_vals = zeros(numel(N_list), 1);
    ratio_vals = nan(numel(N_list), 1);
    nsteps_vals = zeros(numel(N_list), 1);

    for i = 1:numel(N_list)
        N = N_list(i);
        result = burgers_solve(N, nu, T, grid_name, form_name, CFL, []);
        s_vals(i) = result.s_star;
        err_vals(i) = abs(result.s_star - reference_s_star) / reference_s_star;
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
    studies.(key).N = N_list;
    studies.(key).nsteps = nsteps_vals;
    studies.(key).s_star = s_vals;
    studies.(key).rel_err = err_vals;
    studies.(key).ratio = ratio_vals;
    studies.(key).rows = rows;
    figure(fig);
    loglog(N_list, err_vals, ...
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
xlim([min(N_list), max(N_list)]);
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
studies.figure_path = fullfile(root, 'q1e_convergence.png');
print(fig, studies.figure_path, '-dpng', '-r260');
close(fig);
print_summary(studies, case_defs);
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end

function setup_graphics()
setenv('GNUTERM', 'pngcairo');
set(0, 'defaultfigurevisible', 'off');
end

function print_summary(studies, case_defs)
fprintf('\nQ1e results\n');
for c = 1:size(case_defs, 1)
    grid_name = case_defs{c, 1};
    form_name = case_defs{c, 2};
    key = sprintf('%s_%s', grid_name, form_name);
    fprintf('\n%s / %s\n', capitalize_word(grid_name), capitalize_word(form_name));
    fprintf('N\t nsteps\t s*\t\t rel.err\t ratio\n');
    rows = studies.(key).rows;
    for i = 1:size(rows, 1)
        fprintf('%s\t %s\t %s\t %s\t %s\n', rows{i, 1}, rows{i, 2}, rows{i, 3}, rows{i, 4}, rows{i, 5});
    end
end
fprintf('\nPlot saved to %s\n', studies.figure_path);
end
