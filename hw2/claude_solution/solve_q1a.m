clear; clc; close all;

% ==========================================
% Driver for 1D Burgers equation
% ==========================================

% -----------------------------
% Parameters
% -----------------------------
N         = 200;
nu        = 0.01;                 % change if needed
T         = 2.0;
grid_type = 'uniform';            % 'uniform' or 'chebyshev'
form      = 'convective';         % 'convective' or 'conservation'
t_save    = 0:0.1:2.0;

% CFL choice
% For Chebyshev grid, use smaller CFL
if strcmpi(grid_type, 'uniform')
    CFL = 0.1;
else
    CFL = 0.01;
end

% -----------------------------
% Solve
% -----------------------------
results = burgers_solve(N, nu, T, grid_type, form, CFL, t_save);

% Extract results
xb     = results.xb;
u_snap = results.u_snap;
t_snap = results.t_snap;

% -----------------------------
% Plot snapshots on one figure
% -----------------------------
figure;
hold on;

cmap = lines(length(t_snap));
legtxt = cell(length(t_snap), 1);

for k = 1:length(t_snap)
    plot(xb, u_snap(:,k), 'LineWidth', 1.2, 'Color', cmap(k,:));
    legtxt{k} = sprintf('t = %.3f', t_snap(k));
end

xlabel('x', 'FontSize', 12);
ylabel('u(x,t)', 'FontSize', 12);
title(sprintf('1D Burgers Solution, N = %d, %s grid, %s form', ...
      N, grid_type, form), 'FontSize', 12);
xlim([0 1]);
grid on;
legend(legtxt, 'Location', 'eastoutside');
hold off;

% -----------------------------
% Print diagnostics
% -----------------------------
fprintf('Done.\n');
fprintf('Grid type : %s\n', grid_type);
fprintf('Form      : %s\n', form);
fprintf('N         : %d\n', N);
fprintf('dt        : %.6e\n', results.dt);
fprintf('nsteps    : %d\n', results.nsteps);
fprintf('s*        : %.8f\n', results.s_star);
fprintf('t*        : %.8f\n', results.t_star);

% -----------------------------
% Optional: save figure
% -----------------------------
print('-dpng', 'burgers_snapshots.png');
