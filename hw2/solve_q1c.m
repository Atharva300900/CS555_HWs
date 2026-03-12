function study = solve_q1c()
% HW2 Q1c: convergence of s* for the convective form.
setup_graphics();

nu = 1 / (100 * pi);
T = 2.0;
CFL = 0.1;
reference_s_star = 152.00516;
N_list = 2 .^ (5:11);
grids = {'uniform', 'chebyshev'};
root = fileparts(mfilename('fullpath'));
study = struct();

fig = figure('visible', 'off', 'Position', [100 100 900 550]);
hold on;

for g = 1:numel(grids)
    grid_name = grids{g};
    rows = cell(numel(N_list), 5);
    s_vals = zeros(numel(N_list), 1);
    err_vals = zeros(numel(N_list), 1);

    for i = 1:numel(N_list)
        N = N_list(i);
        result = burgers_solve(N, nu, T, grid_name, 'convective', CFL, []);
        s_vals(i) = result.s_star;
        err_vals(i) = abs(result.s_star - reference_s_star) / reference_s_star;
        rows(i, :) = {sprintf('%d', N), sprintf('%.6f', result.dt), sprintf('%d', result.nsteps), ...
            sprintf('%.8f', result.s_star), sprintf('%.3e', err_vals(i))};
    end

    study.(grid_name).N = N_list;
    study.(grid_name).s_star = s_vals;
    study.(grid_name).rel_err = err_vals;
    study.(grid_name).five_digit_estimate = sprintf('%.5g', s_vals(end));
    study.(grid_name).rows = rows;

    figure(fig);
    loglog(N_list, err_vals, 'o-', 'LineWidth', 1.4, 'MarkerSize', 7);
end

figure(fig);
hold off;
grid on;
xlabel('N');
ylabel('relative error in s^*');
title('Q1c: convergence of s^* for the convective form');
legend({'Uniform', 'Chebyshev'}, 'Location', 'southwest');
study.figure_path = fullfile(root, 'q1c_convergence.png');
print(fig, study.figure_path, '-dpng', '-r200');
close(fig);

study.reference_s_star = reference_s_star;
study.reference_s_star_5sf = sprintf('%.5g', reference_s_star);
print_summary(study, grids);
end

function out = capitalize_word(str)
out = lower(str);
out(1) = upper(out(1));
end

function setup_graphics()
setenv('GNUTERM', 'pngcairo');
set(0, 'defaultfigurevisible', 'off');
end

function print_summary(study, grids)
fprintf('\nQ1c results\n');
fprintf('Reference s* = %.8f\n', study.reference_s_star);
for g = 1:numel(grids)
    grid_name = grids{g};
    fprintf('\n%s grid\n', capitalize_word(grid_name));
    fprintf('N\t dt\t\t nsteps\t s*\t\t rel.err\n');
    rows = study.(grid_name).rows;
    for i = 1:size(rows, 1)
        fprintf('%s\t %s\t %s\t %s\t %s\n', rows{i, 1}, rows{i, 2}, rows{i, 3}, rows{i, 4}, rows{i, 5});
    end
    fprintf('5 significant digits: %s\n', study.(grid_name).five_digit_estimate);
end
fprintf('\nPlot saved to %s\n', study.figure_path);
end
