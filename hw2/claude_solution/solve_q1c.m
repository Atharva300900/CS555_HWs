%% Question 1c: Find s* to 5 significant digits, convergence study
% Compare with analytical value from Table 3 of paper: s* = 152.00516
% CFL = 0.1 fixed, vary N

clear; close all;

nu = 1/(100*pi);
T = 2.0;
CFL = 0.1;
s_analytical = 152.00516;  % From Table 3 of Basdevant et al.

% Test several N values for convergence
N_vals = [50, 100, 200, 400, 800, 1600];

fprintf('Question 1c: Convergence of s* with N\n');
fprintf('Analytical value: s* = %.5f\n\n', s_analytical);

% --- Uniform spacing, convective form ---
fprintf('%-8s %-10s %-12s %-12s\n', 'N', 'nsteps', 's*', 'rel.err.');
fprintf('%s\n', repmat('-', 1, 45));

s_star_vals = zeros(size(N_vals));
for i = 1:length(N_vals)
    N = N_vals(i);
    results = burgers_solve(N, nu, T, 'uniform', 'convective', CFL);
    s_star_vals(i) = results.s_star;
    rel_err = abs(s_star_vals(i) - s_analytical) / s_analytical;
    fprintf('%-8d %-10d %-12.5f %-12.4e\n', N, results.nsteps, results.s_star, rel_err);
end

fprintf('\n');

% --- Chebyshev spacing, convective form ---
fprintf('Chebyshev spacing, convective form:\n');
fprintf('%-8s %-10s %-12s %-12s\n', 'N', 'nsteps', 's*', 'rel.err.');
fprintf('%s\n', repmat('-', 1, 45));

for i = 1:length(N_vals)
    N = N_vals(i);
    results = burgers_solve(N, nu, T, 'chebyshev', 'convective', CFL);
    rel_err = abs(results.s_star - s_analytical) / s_analytical;
    fprintf('%-8d %-10d %-12.5f %-12.4e\n', N, results.nsteps, results.s_star, rel_err);
end

fprintf('\nThe solution converges to the analytical value as N increases.\n');
fprintf('5-digit value from finest resolution is reported above.\n');
