%% Question 1e: Convergence tables for four cases
% N = 2^k, k = 5, 6, ..., 11
% Four cases: {uniform, chebyshev} x {convective, conservation}

clear; close all;

nu = 1/(100*pi);
T = 2.0;
CFL = 0.1;
s_analytical = 152.00516;  % From Table 3 of Basdevant et al.

k_vals = 5:11;
N_vals = 2.^k_vals;  % [32, 64, 128, 256, 512, 1024, 2048]

cases = {
    'uniform',   'convective',   'Uniform spacing, convective form';
    'chebyshev', 'convective',   'Chebyshev spacing, convective form';
    'uniform',   'conservation', 'Uniform spacing, conservation form';
    'chebyshev', 'conservation', 'Chebyshev spacing, conservation form';
};

% Store all results for latex table generation
all_results = cell(size(cases, 1), 1);

for c = 1:size(cases, 1)
    grid_type = cases{c, 1};
    form_type = cases{c, 2};
    case_name = cases{c, 3};

    fprintf('\n========================================\n');
    fprintf('Case %d: %s\n', c, case_name);
    fprintf('========================================\n');
    fprintf('%-8s %-10s %-12s %-12s %-8s\n', 'N', 'nsteps', 's*', 'rel.err.', 'ratio');
    fprintf('%s\n', repmat('-', 1, 55));

    s_star_vals = zeros(size(N_vals));
    err_vals = zeros(size(N_vals));
    nsteps_vals = zeros(size(N_vals));

    for i = 1:length(N_vals)
        N = N_vals(i);
        results = burgers_solve(N, nu, T, grid_type, form_type, CFL);
        s_star_vals(i) = results.s_star;
        err_vals(i) = abs(results.s_star - s_analytical) / s_analytical;
        nsteps_vals(i) = results.nsteps;

        if i == 1
            ratio_str = '---';
        else
            ratio = err_vals(i-1) / err_vals(i);
            ratio_str = sprintf('%.2f', ratio);
        end

        fprintf('%-8d %-10d %-12.5f %-12.4e %-8s\n', ...
            N, results.nsteps, results.s_star, err_vals(i), ratio_str);
    end

    all_results{c}.N = N_vals;
    all_results{c}.nsteps = nsteps_vals;
    all_results{c}.s_star = s_star_vals;
    all_results{c}.err = err_vals;
    all_results{c}.name = case_name;
end

% Save results to .mat file for LaTeX table generation
save('q1e_results.mat', 'all_results', 'N_vals', 's_analytical');

fprintf('\n\nDiscussion:\n');
fprintf('- For uniform spacing, the convergence ratio should approach 4\n');
fprintf('  (2nd-order convergence in space) as N increases.\n');
fprintf('- Chebyshev spacing clusters points near boundaries x=0 and x=1,\n');
fprintf('  which helps resolve the steep gradient near x=1, potentially\n');
fprintf('  giving better accuracy for the same N.\n');
fprintf('- The convective and conservation forms are formally equivalent\n');
fprintf('  but their discrete approximations differ. The convective form\n');
fprintf('  typically gives slightly better results for this problem.\n');
fprintf('- At coarse N, the ratio may deviate from 4 due to under-resolution\n');
fprintf('  of the thin internal layer.\n');

% Generate convergence plot
figure('Position', [100 100 800 500]);
markers = {'o-', 's-', 'o--', 's--'};
colors = {'b', 'r', 'b', 'r'};
for c = 1:4
    loglog(all_results{c}.N, all_results{c}.err, ...
        markers{c}, 'Color', colors{c}, 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
end
% Reference line for 2nd order
loglog(N_vals, 2*N_vals(1)^2 ./ N_vals.^2 * all_results{1}.err(1), ...
    'k--', 'LineWidth', 1, 'DisplayName', 'O(N^{-2})');

xlabel('N', 'FontSize', 14);
ylabel('Relative Error in s^*', 'FontSize', 14);
title('Convergence of s^* for Four Cases', 'FontSize', 14);
legend(cases{1,3}, cases{2,3}, cases{3,3}, cases{4,3}, 'O(N^{-2})', ...
    'Location', 'southwest', 'FontSize', 10);
grid on;
saveas(gcf, 'q1e_convergence.png');
fprintf('\nFigure saved to q1e_convergence.png\n');
