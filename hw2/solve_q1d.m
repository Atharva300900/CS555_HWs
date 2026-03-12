function summary = solve_q1d()
% HW2 Q1d: compare convective and conservation forms at the finest N.
setup_graphics();

nu = 1 / (100 * pi);
T = 2.0;
CFL = 0.1;
reference_s_star = 152.00516;
N_list = 2 .^ (5:11);
N = N_list(end);
root = fileparts(mfilename('fullpath'));
grids = {'uniform', 'chebyshev'};

summary = struct();
fig = figure('visible', 'off', 'Position', [100 100 950 650]);

for i = 1:numel(grids)
    grid_name = grids{i};
    convective = burgers_solve(N, nu, T, grid_name, 'convective', CFL, []);
    conservative = burgers_solve(N, nu, T, grid_name, 'conservation', CFL, []);

    summary.(grid_name).convective = convective;
    summary.(grid_name).conservative = conservative;
    summary.(grid_name).difference = abs(convective.s_star - conservative.s_star);
    summary.(grid_name).rel_err_conservation = abs(conservative.s_star - reference_s_star) / reference_s_star;

    subplot(2, 1, i);
    plot(convective.t_all, convective.s_all, 'b-', 'LineWidth', 1.3); hold on;
    plot(conservative.t_all, conservative.s_all, 'r--', 'LineWidth', 1.3);
    plot([0.0, T], [reference_s_star, reference_s_star], 'k:', 'LineWidth', 1.0);
    hold off;
    grid on;
    xlabel('t');
    ylabel('s(t)');
    title(sprintf('Q1d: %s grid, N = %d', capitalize_word(grid_name), N));
    legend({'Convective', 'Conservation', 'Reference s^*'}, 'Location', 'southeast');
end

summary.figure_path = fullfile(root, 'q1d_form_comparison.png');
print(fig, summary.figure_path, '-dpng', '-r200');
close(fig);

print_summary(summary, grids);
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end

function setup_graphics()
setenv('GNUTERM', 'pngcairo');
set(0, 'defaultfigurevisible', 'off');
end

function print_summary(summary, grids)
fprintf('\nQ1d results\n');
for i = 1:numel(grids)
    grid_name = grids{i};
    convective = summary.(grid_name).convective;
    conservative = summary.(grid_name).conservative;
    fprintf('%-10s convective=%.8f conservation=%.8f rel.err(cons)=%.3e diff=%.3e\n', ...
        capitalize_word(grid_name), convective.s_star, conservative.s_star, ...
        summary.(grid_name).rel_err_conservation, summary.(grid_name).difference);
end
fprintf('Plot saved to %s\n', summary.figure_path);
end
