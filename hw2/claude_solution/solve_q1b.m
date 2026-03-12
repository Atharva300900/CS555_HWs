clear; clc; close all;

% ==========================================
% Q1(b): slope history s(t) = max_x |du/dx|
% ==========================================

% -----------------------------
% User settings
% -----------------------------
N         = 200;
nu        = 0.01;              % change if needed
T         = 2.0;
grid_type = 'chebyshev';         % 'uniform' or 'chebyshev'
form      = 'convective';      % 'convective' or 'conservation'

% Baseline CFL
if strcmpi(grid_type,'uniform')
    CFL = 0.10;
else
    CFL = 0.01;
end

% Optional reference values from analytical/benchmark result
use_reference = true;
s_ref = 152.00516;
t_ref = 0.5055;

% -----------------------------
% Baseline run
% -----------------------------
results = burgers_solve(N, nu, T, grid_type, form, CFL);

% -----------------------------
% Verification 1: dt/2
% -----------------------------
results_dt2 = burgers_solve(N, nu, T, grid_type, form, CFL/2);

% -----------------------------
% Verification 2: finer grid
% -----------------------------
N2 = 2*N;

% For uniform grid, increase CFL so dt stays about the same:
% dt = CFL/N, so for 2N use CFL2 = 2*CFL
if strcmpi(grid_type,'uniform')
    CFL2 = 2*CFL;
else
    CFL2 = CFL;   % for Chebyshev just keep same CFL
end

results_N2 = burgers_solve(N2, nu, T, grid_type, form, CFL2);

% -----------------------------
% Print diagnostics
% -----------------------------
fprintf('\n========================================\n');
fprintf('Q1(b): slope history s(t) = max_x |du/dx|\n');
fprintf('========================================\n');
fprintf('Grid type         : %s\n', grid_type);
fprintf('Form              : %s\n', form);
fprintf('N                  = %d\n', N);
fprintf('nu                 = %.8g\n', nu);
fprintf('T                  = %.8g\n', T);
fprintf('Baseline CFL       = %.8g\n', CFL);
fprintf('Baseline dt        = %.8e\n', results.dt);
fprintf('Baseline nsteps    = %d\n', results.nsteps);
fprintf('\n');
fprintf('Baseline result:\n');
fprintf('s*                 = %.10f\n', results.s_star);
fprintf('t*                 = %.10f\n', results.t_star);
fprintf('\n');
fprintf('Verification with dt/2:\n');
fprintf('dt(dt/2 run)       = %.8e\n', results_dt2.dt);
fprintf('s*(dt/2)           = %.10f\n', results_dt2.s_star);
fprintf('t*(dt/2)           = %.10f\n', results_dt2.t_star);
fprintf('relative change s* = %.6e\n', ...
    abs(results_dt2.s_star - results.s_star) / max(abs(results_dt2.s_star),1e-14));
fprintf('absolute change t* = %.6e\n', ...
    abs(results_dt2.t_star - results.t_star));
fprintf('\n');
fprintf('Verification with finer grid:\n');
fprintf('N_fine             = %d\n', N2);
fprintf('dt(fine run)       = %.8e\n', results_N2.dt);
fprintf('s*(fine grid)      = %.10f\n', results_N2.s_star);
fprintf('t*(fine grid)      = %.10f\n', results_N2.t_star);
fprintf('relative change s* = %.6e\n', ...
    abs(results_N2.s_star - results.s_star) / max(abs(results_N2.s_star),1e-14));
fprintf('absolute change t* = %.6e\n', ...
    abs(results_N2.t_star - results.t_star));
fprintf('========================================\n\n');

% -----------------------------
% Plot baseline and dt/2 check
% -----------------------------
figure;
hold on;
plot(results.t_all, results.s_all, 'b-', 'LineWidth', 2);
plot(results_dt2.t_all, results_dt2.s_all, 'r--', 'LineWidth', 1.3);

% Mark numerical maximum
plot(results.t_star, results.s_star, 'ro', 'MarkerSize', 12, 'LineWidth', 2);

% Optional reference lines
if use_reference
    yline(s_ref, 'k--', 'LineWidth', 1.2);
    xline(t_ref, 'k:',  'LineWidth', 1.5);
end

xlabel('t', 'FontSize', 14);
ylabel('s(t) = max |u_x|', 'FontSize', 14);
title(sprintf('Q1(b): slope history for N = %d, %s form, %s grid', ...
    N, form, grid_type), 'FontSize', 16);
grid on;
xlim([0 T]);

% Annotation near the peak
text(results.t_star + 0.05, results.s_star - 8, ...
    sprintf('s^* = %.4f\\nt^* = %.4f', results.s_star, results.t_star), ...
    'FontSize', 12);

if use_reference
    legend('baseline', 'dt/2 check', 'numerical peak', 'reference s^*', 'reference t^*', ...
        'Location', 'best');
else
    legend('baseline', 'dt/2 check', 'numerical peak', 'Location', 'best');
end
hold off;

print('-dpng', 'q1b_slope_history.png');

% -----------------------------
% Separate plot: fine-grid comparison
% -----------------------------
figure;
hold on;
plot(results.t_all, results.s_all, 'b-', 'LineWidth', 2);
plot(results_N2.t_all, results_N2.s_all, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8);

if use_reference
    yline(s_ref, 'k--', 'LineWidth', 1.2);
    xline(t_ref, 'k:',  'LineWidth', 1.5);
end

xlabel('t', 'FontSize', 14);
ylabel('s(t) = max |u_x|', 'FontSize', 14);
title(sprintf('Q1(b): grid verification, %s form, %s grid', form, grid_type), ...
    'FontSize', 16);
grid on;
xlim([0 T]);

if use_reference
    legend(sprintf('N = %d',N), sprintf('N = %d',N2), 'reference s^*', 'reference t^*', ...
        'Location', 'best');
else
    legend(sprintf('N = %d',N), sprintf('N = %d',N2), 'Location', 'best');
end
hold off;

print('-dpng', 'q1b_grid_verification.png');
