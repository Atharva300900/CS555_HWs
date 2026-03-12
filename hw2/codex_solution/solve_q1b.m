function summary = solve_q1b()
% HW2 Q1b: compute s(t), locate its maximum, and verify with dt/2.
hw2_setup_graphics();
cfg = hw2_config();
paths = hw2_paths();

N = cfg.N_q1a;
target_dt = cfg.fixed_cfl / N;
grids = {'uniform', 'chebyshev'};
summary = struct();
rows = cell(numel(grids), 8);

fig_t = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;
fig_pit = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;

for i = 1:numel(grids)
    grid_name = grids{i};
    result = hw2_load_or_solve(N, grid_name, 'convective', target_dt, []);
    verify = hw2_load_or_solve(N, grid_name, 'convective', target_dt / 2.0, []);
    rel_change = abs(result.s_star - verify.s_star) / max(abs(verify.s_star), eps);
    rel_err = abs(result.s_star - cfg.reference_s_star) / cfg.reference_s_star;

    summary.(grid_name).base = result;
    summary.(grid_name).half_dt = verify;
    summary.(grid_name).rel_change_dt_half = rel_change;
    summary.(grid_name).rel_err_reference = rel_err;

    figure(fig_t);
    plot(result.t_all, result.s_all, 'LineWidth', 1.4);
    figure(fig_pit);
    plot(pi * result.t_all, result.s_all, 'LineWidth', 1.4);

    rows(i, :) = {capitalize_word(grid_name), sprintf('%.6f', result.dt), sprintf('%d', result.nsteps), ...
        sprintf('%.8f', result.s_star), sprintf('%.8f', result.t_star), sprintf('%.3e', rel_err), ...
        sprintf('%.6f', verify.dt), sprintf('%.3e', rel_change)};
end

figure(fig_t);
y_lim = ylim;
plot([0.0, cfg.T], [cfg.reference_s_star, cfg.reference_s_star], 'k--', 'LineWidth', 1.0);
plot([cfg.reference_t_star, cfg.reference_t_star], y_lim, 'k:', 'LineWidth', 1.0);
ylim(y_lim);
grid on;
xlabel('t');
ylabel('s(t) = max |u_x|');
title('Q1b: s(t) versus t for N = 200');
legend({'Uniform', 'Chebyshev', 'Reference s^*', 'Reference t^*'}, 'Location', 'southeast');
summary.figure_t = fullfile(paths.figures, 'q1b_s_vs_t.png');
print(fig_t, summary.figure_t, '-dpng', '-r200');
close(fig_t);

figure(fig_pit);
y_lim = ylim;
plot([0.0, pi * cfg.T], [cfg.reference_s_star, cfg.reference_s_star], 'k--', 'LineWidth', 1.0);
plot([cfg.reference_t_star_paper, cfg.reference_t_star_paper], y_lim, 'k:', 'LineWidth', 1.0);
ylim(y_lim);
grid on;
xlabel('pi t');
ylabel('s(t) = max |u_x|');
title('Q1b: s(t) versus pi t for comparison with Fig. 6');
legend({'Uniform', 'Chebyshev', 'Reference s^*', 'Reference pi t^*'}, 'Location', 'southeast');
summary.figure_pit = fullfile(paths.figures, 'q1b_s_vs_pit.png');
print(fig_pit, summary.figure_pit, '-dpng', '-r200');
close(fig_pit);

summary.table_path = fullfile(paths.tables, 'q1b_summary.txt');
hw2_write_table_txt(summary.table_path, ...
    {'Grid', 'dt', 'nsteps', 's_star', 't_star', 'rel_err_ref', 'dt_half', 'rel_change_dt_half'}, rows);
save(fullfile(paths.data, 'q1b_summary.mat'), 'summary');
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end
