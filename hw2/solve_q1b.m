function summary = solve_q1b()
% HW2 Q1b: compute s(t), locate its maximum, and verify with dt/2.
setup_graphics();

nu = 1 / (100 * pi);
T = 2.0;
N = 200;
CFL = 0.1;
reference_s_star = 152.00516;
reference_t_star_paper = 1.6037;
reference_t_star = reference_t_star_paper / pi;
root = fileparts(mfilename('fullpath'));
grids = {'uniform', 'chebyshev'};

summary = struct();
plot_data = struct('t', [], 'pi_t', [], 's_uniform', [], 's_chebyshev', []);

fig_t = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;

for i = 1:numel(grids)
    grid_name = grids{i};
    result = burgers_solve(N, nu, T, grid_name, 'convective', CFL, []);
    verify = burgers_solve(N, nu, T, grid_name, 'convective', CFL / 2.0, []);
    rel_change = abs(result.s_star - verify.s_star) / max(abs(verify.s_star), eps);
    rel_err = abs(result.s_star - reference_s_star) / reference_s_star;

    summary.(grid_name).base = result;
    summary.(grid_name).half_dt = verify;
    summary.(grid_name).rel_change_dt_half = rel_change;
    summary.(grid_name).rel_err_reference = rel_err;

    if isempty(plot_data.t)
        plot_data.t = result.t_all(:);
        plot_data.pi_t = pi * result.t_all(:);
    end
    plot_data.(sprintf('s_%s', lower(grid_name))) = result.s_all(:);

    figure(fig_t);
    plot(result.t_all, result.s_all, 'LineWidth', 1.4);
end

figure(fig_t);
y_lim = ylim;
plot([0.0, T], [reference_s_star, reference_s_star], 'k--', 'LineWidth', 1.0);
plot([reference_t_star, reference_t_star], y_lim, 'k:', 'LineWidth', 1.0);
ylim(y_lim);
grid on;
xlabel('t');
ylabel('s(t) = max |u_x|');
title('Q1b: s(t) versus t for N = 200');
legend({'Uniform', 'Chebyshev', 'Reference s^*', 'Reference t^*'}, 'Location', 'southeast');
summary.figure_t = fullfile(root, 'q1b_s_vs_t.png');
print(fig_t, summary.figure_t, '-dpng', '-r200');
close(fig_t);

summary.plot_data = plot_data;
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
fprintf('\nQ1b results\n');
for i = 1:numel(grids)
    grid_name = grids{i};
    result = summary.(grid_name).base;
    fprintf('%-10s dt=%0.6f nsteps=%d s*=%.8f t*=%.8f rel.err=%0.3e dt/2 change=%0.3e\n', ...
        capitalize_word(grid_name), result.dt, result.nsteps, result.s_star, result.t_star, ...
        summary.(grid_name).rel_err_reference, summary.(grid_name).rel_change_dt_half);
end
fprintf('Plot saved to %s\n', summary.figure_t);
fprintf('Raw plot values are returned in summary.plot_data.\n');
end
