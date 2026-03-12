function study = solve_q1c()
% HW2 Q1c: convergence of s* for the convective form.
hw2_setup_graphics();
cfg = hw2_config();
paths = hw2_paths();
grids = {'uniform', 'chebyshev'};
study = struct();

fig = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;
fig_ref = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;

for g = 1:numel(grids)
    grid_name = grids{g};
    rows = cell(numel(cfg.N_list), 5);
    s_vals = zeros(numel(cfg.N_list), 1);
    err_vals = zeros(numel(cfg.N_list), 1);

    for i = 1:numel(cfg.N_list)
        N = cfg.N_list(i);
        target_dt = cfg.fixed_cfl / N;
        result = hw2_load_or_solve(N, grid_name, 'convective', target_dt, []);
        s_vals(i) = result.s_star;
        err_vals(i) = abs(result.s_star - cfg.reference_s_star) / cfg.reference_s_star;
        rows(i, :) = {sprintf('%d', N), sprintf('%.6f', result.dt), sprintf('%d', result.nsteps), ...
            sprintf('%.8f', result.s_star), sprintf('%.3e', err_vals(i))};
    end

    study.(grid_name).N = cfg.N_list;
    study.(grid_name).s_star = s_vals;
    study.(grid_name).rel_err = err_vals;
    study.(grid_name).five_digit_estimate = sprintf('%.5g', s_vals(end));

    hw2_write_table_txt(fullfile(paths.tables, sprintf('q1c_%s.txt', grid_name)), ...
        {'N', 'dt', 'nsteps', 's_star', 'rel_err'}, rows);
    figure(fig);
    loglog(cfg.N_list, err_vals, 'o-', 'LineWidth', 1.4, 'MarkerSize', 7);
    figure(fig_ref);
    semilogx(cfg.N_list, s_vals, 'o-', 'LineWidth', 1.4, 'MarkerSize', 7);
end

figure(fig);
hold off;
grid on;
xlabel('N');
ylabel('relative error in s^*');
title('Q1c: convergence of s^* for the convective form');
legend({'Uniform', 'Chebyshev'}, 'Location', 'southwest');
study.figure_path = fullfile(paths.figures, 'q1c_convergence.png');
print(fig, study.figure_path, '-dpng', '-r200');
close(fig);

figure(fig_ref);
plot(cfg.N_list, cfg.reference_s_star * ones(size(cfg.N_list)), 'k--', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('N');
ylabel('s^*');
title('Q1c: computed s^* versus analytical reference');
legend({'Uniform', 'Chebyshev', 'Analytical s^*'}, 'Location', 'southeast');
study.reference_figure_path = fullfile(paths.figures, 'q1c_s_star_vs_reference.png');
print(fig_ref, study.reference_figure_path, '-dpng', '-r200');
close(fig_ref);

study.reference_s_star = cfg.reference_s_star;
study.reference_s_star_5sf = sprintf('%.5g', cfg.reference_s_star);
save(fullfile(paths.data, 'q1c_study.mat'), 'study');
end
