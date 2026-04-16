function hw3()
% CS555 HW4: Advection-Diffusion with Linear FEM + BDF3/EXT3
%
%   u_t + c u_x = nu u_xx + f   on [0, L]
%
% Q1: Steady state (u_t=0), f=1, u(0)=u(1)=0
  % Q2: Unsteady,     f=0,    u(0,t)=sin(pi*t)

  if any(strcmp(available_graphics_toolkits(), 'gnuplot'))
    graphics_toolkit('gnuplot');
  end
  set(0, 'defaultfigurevisible', 'off');

  L = 1;  c = 1;  nu = 1e-3;  f = 1;

  fprintf('============================================\n');
  fprintf('  CS555 HW4 — Advection-Diffusion FEM\n');
  fprintf('============================================\n\n');

  %% =============== PROBLEM 1: STEADY STATE ===============
  fprintf('---------- Problem 1: Steady State ----------\n\n');

  % --- Q1(a): uniform mesh, E = 100 ---
  E = 100;
  x = linspace(0, L, E+1)';
  [u_a, ue_a] = solve_steady(x, c, nu, f, L);
  [abs_a, rel_a] = compute_error(u_a, ue_a);
  fprintf('Q1(a): E = %d, uniform mesh\n', E);
  fprintf('  Max pointwise relative error: %.6e\n', rel_a);
  fprintf('  Max absolute error:           %.6e\n\n', abs_a);
  plot_q1(x, u_a, c, nu, f, L, abs_a, ...
          'Q1(a): uniform, E = 100', 'q1a_uniform.png');

  % --- Q1(b): geometric mesh, E = 20, s = 0.7 ---
  E = 20;  s = 0.7;
  x = geometric_mesh(L, E, s);
  [u_b, ue_b] = solve_steady(x, c, nu, f, L);
  [abs_b, rel_b] = compute_error(u_b, ue_b);
  fprintf('Q1(b): E = %d, s = %.2f, geometric mesh\n', E, s);
  fprintf('  Max pointwise relative error: %.6e\n', rel_b);
  fprintf('  Max absolute error:           %.6e\n\n', abs_b);
  plot_q1(x, u_b, c, nu, f, L, abs_b, ...
          sprintf('Q1(b): geometric, E = %d, s = %.2f', E, s), ...
          'q1b_geometric_s0p7.png');

  % --- Q1(c): find optimal s for E = 20 ---
  E = 20;
  s_vals = 0.20:0.001:1.05;
  abs_vs_s = zeros(size(s_vals));
  rel_vs_s = zeros(size(s_vals));
  for k = 1:length(s_vals)
    xk = geometric_mesh(L, E, s_vals(k));
    [uk, uek] = solve_steady(xk, c, nu, f, L);
    [abs_vs_s(k), rel_vs_s(k)] = compute_error(uk, uek);
  end
  [best_rel, idx] = min(rel_vs_s);
  s_opt = s_vals(idx);

  x = geometric_mesh(L, E, s_opt);
  [u_c, ue_c] = solve_steady(x, c, nu, f, L);
  [abs_c, rel_c] = compute_error(u_c, ue_c);
  fprintf('Q1(c): E = %d, optimal s search in [%.2f, %.2f]\n', ...
          E, s_vals(1), s_vals(end));
  fprintf('  Metric: pointwise relative error\n');
  fprintf('  Best s:                       %.3f\n', s_opt);
  fprintf('  Max pointwise relative error: %.6e\n', rel_c);
  fprintf('  Max absolute error:           %.6e\n\n', abs_c);
  plot_q1(x, u_c, c, nu, f, L, abs_c, ...
          sprintf('Q1(c): optimal, E = %d, s = %.3f', E, s_opt), ...
          'q1c_best_s.png');

  % Error-vs-s plot
  fig = figure('visible', 'off');
  semilogy(s_vals, rel_vs_s, 'b-o', 'linewidth', 1.5, 'markersize', 3);
  hold on;
  semilogy(s_opt, best_rel, 'rp', 'markersize', 14, 'linewidth', 2);
  hold off;
  xlabel('s');  ylabel('Max pointwise relative error');
  title('Q1(c): Relative error vs geometric ratio s  (E = 20)');
  legend('Relative error(s)', sprintf('Best s = %.3f', s_opt), ...
         'location', 'northeast');
  grid on;
  print(fig, 'q1c_error_vs_s.png', '-dpng', '-r200');
  close(fig);

  %% =============== PROBLEM 2: UNSTEADY ===============
  fprintf('---------- Problem 2: Unsteady BDF3/EXT3 ----------\n\n');

  E  = 100;
  x  = linspace(0, L, E+1)';
  n  = E + 1;
  h  = L / E;
  dt = 1e-3;                   % CFL = c*dt/h = 0.1
  t_end  = 1.5;
  nsteps = round(t_end / dt);
  dt = t_end / nsteps;        % adjust so we hit t_end exactly

  [Bbar, Cbar, Abar] = assemble_fem(x, c);

  % --- Q2(a): Dirichlet both ends ---
  ua = bdf3ext3(Bbar, Cbar, Abar, n, dt, nsteps, nu, ...
                @(t) sin(pi*t), true, @(t) 0);
  fprintf('Q2(a): Dirichlet BCs,  E = %d,  dt = %.4e,  %d steps\n', E, dt, nsteps);
  fprintf('  min(u) = %+.6e,  max(u) = %+.6e,  u(1) = %+.6e\n\n', ...
          min(ua), max(ua), ua(end));

  fig = figure('visible', 'off');
  plot(x, ua, 'b-', 'linewidth', 2);
  xlabel('x');  ylabel('u(x, t)');
  title(sprintf('Q2(a): t = %.2f,  u(1, t) = 0', t_end));
  grid on;
  print(fig, 'q2a_t1p5.png', '-dpng', '-r200');
  close(fig);

  % --- Q2(b): Left Dirichlet, Right Neumann ---
  ub = bdf3ext3(Bbar, Cbar, Abar, n, dt, nsteps, nu, ...
                @(t) sin(pi*t), false, []);
  fprintf('Q2(b): Neumann at x=L, E = %d,  dt = %.4e,  %d steps\n', E, dt, nsteps);
  fprintf('  min(u) = %+.6e,  max(u) = %+.6e,  u(1) = %+.6e\n\n', ...
          min(ub), max(ub), ub(end));

  fig = figure('visible', 'off');
  plot(x, ub, 'b-', 'linewidth', 2);
  xlabel('x');  ylabel('u(x, t)');
  title(sprintf('Q2(b): t = %.2f,  du/dx(1, t) = 0', t_end));
  grid on;
  print(fig, 'q2b_t1p5.png', '-dpng', '-r200');
  close(fig);

  % Comparison plot
  fig = figure('visible', 'off');
  plot(x, ua, 'b-', 'linewidth', 2);  hold on;
  plot(x, ub, 'r--', 'linewidth', 2); hold off;
  xlabel('x');  ylabel('u(x, t)');
  title(sprintf('Q2: Comparison at t = %.2f', t_end));
  legend('(a) u(1,t) = 0', '(b) du/dx(1,t) = 0', 'location', 'northwest');
  grid on;
  print(fig, 'q2_comparison.png', '-dpng', '-r200');
  close(fig);

  %% =============== WRITE SUMMARY FILE ===============
  write_results(abs_a, rel_a, abs_b, rel_b, s, abs_c, rel_c, s_opt, ...
                ua, ub, x, c, nu, t_end);
  fprintf('All plots saved.  Summary in hw3_results.txt.\n');
end


%% =====================================================================
%%  ASSEMBLY:  Abar = Q' * A_L * Q   (and similarly B, C)
%% =====================================================================
function [Bbar, Cbar, Abar] = assemble_fem(xnodes, c)
% 1-D linear FEM assembly via the scatter operator Q.
%
% Reference element  r in [-1, 1]  with local bases
%   l_0(r) = (1-r)/2,    l_1(r) = (1+r)/2
%   dl_0/dr = -1/2,      dl_1/dr = +1/2
%
% Element mapping: x^e(r) = x_{e-1} + L_e/2 * (r+1),  dx/dr = L_e/2
%
% Local matrices (exact 2-point quadrature, no approximation):
%   A^e = (1/L_e) [1 -1; -1  1]     stiffness
%   B^e = (L_e/6) [2  1;  1  2]     mass
%   C^e = (c/2)   [-1 1; -1  1]     advection  (independent of L_e)

  E  = length(xnodes) - 1;          % number of elements
  n  = E + 1;                        % number of global nodes (0..N)
  nL = 2 * E;                        % total local DOFs

  % ---- scatter operator Q  (nL x n, sparse Boolean) ----
  %  row l = 2*(e-1)+v+1  (v=0,1),   col i = l2g(v,e) = e+v  (1-indexed)
  Qi = zeros(nL, 1);   Qj = zeros(nL, 1);
  for e = 1:E
    Qi(2*e - 1) = 2*e - 1;   Qj(2*e - 1) = e;       % local node 0
    Qi(2*e)     = 2*e;       Qj(2*e)     = e + 1;   % local node 1
  end
  Q = sparse(Qi, Qj, 1, nL, n);

  % ---- block-diagonal local matrices  (nL x nL) ----
  ntriplets = 4 * E;
  ri = zeros(ntriplets, 1);  ci = zeros(ntriplets, 1);
  av = zeros(ntriplets, 1);  bv = zeros(ntriplets, 1);  cv = zeros(ntriplets, 1);

  k = 0;
  for e = 1:E
    Le = xnodes(e + 1) - xnodes(e);

    Ae = (1 / Le) * [1, -1; -1, 1];
    Be = (Le / 6) * [2,  1;  1, 2];
    Ce = (c / 2)  * [-1, 1; -1, 1];

    base = 2 * (e - 1);
    for p = 1:2
      for q = 1:2
        k = k + 1;
        ri(k) = base + p;   ci(k) = base + q;
        av(k) = Ae(p, q);   bv(k) = Be(p, q);   cv(k) = Ce(p, q);
      end
    end
  end

  AL = sparse(ri, ci, av, nL, nL);
  BL = sparse(ri, ci, bv, nL, nL);
  CL = sparse(ri, ci, cv, nL, nL);

  % ---- global matrices:  Q^T * (block-diag) * Q ----
  Abar = Q' * AL * Q;
  Bbar = Q' * BL * Q;
  Cbar = Q' * CL * Q;
end


%% =====================================================================
%%  GEOMETRIC MESH
%% =====================================================================
function x = geometric_mesh(L, E, s)
% L_e = s * L_{e-1},  e = 2..E.   Sum(L_e) = L.
  if abs(s - 1) < 1e-14
    x = linspace(0, L, E + 1)';
  else
    L1 = L * (1 - s) / (1 - s^E);
    Le = L1 * s .^ (0:E-1)';
    x  = [0; cumsum(Le)];
    x(end) = L;                      % enforce exact endpoint
  end
end


%% =====================================================================
%%  EXACT STEADY SOLUTION
%% =====================================================================
function u = exact_steady(x, c, nu, f, L)
% Solves  c u' = nu u'' + f,   u(0) = u(L) = 0.
%
% u(x) = (f/c) x  -  (fL/c) [ exp(c(x-L)/nu) - exp(-cL/nu) ]
%                              / [ 1 - exp(-cL/nu) ]
%
% This form avoids overflow from exp(cL/nu).
  eta = exp(-c * L / nu);                        % ~ 0 when cL/nu >> 1
  u   = (f / c) * x ...
      - (f * L / c) * (exp(c * (x - L) / nu) - eta) ./ (1 - eta);
end


%% =====================================================================
%%  STEADY SOLVER  (restricted system)
%% =====================================================================
function [u, ue] = solve_steady(x, c, nu, f, L)
% Weak form (steady, homogeneous Dirichlet):
%   (nu A + C) u_int = R fbar
%
% where  A = R Abar R',  C = R Cbar R',  fbar_i = int phi_i f dx = f*Bbar*1

  [Bbar, Cbar, Abar] = assemble_fem(x, c);
  n = length(x);

  % Restriction R: keeps interior DOFs  (nodes 2..n-1)
  R = build_restriction(n, true, true);
  Rt = R';

  A = R * Abar * Rt;
  C = R * Cbar * Rt;

  % Forcing vector:  fbar = f * Bbar * ones   (since f is constant)
  fbar = f * Bbar * ones(n, 1);

  u_int = (nu * A + C) \ (R * fbar);

  u = Rt * u_int;

  ue = exact_steady(x, c, nu, f, L);
end


%% =====================================================================
%%  ERROR METRICS
%% =====================================================================
function [max_abs, max_rel] = compute_error(u, ue)
  max_abs = max(abs(u - ue));

  % Pointwise relative error at interior nodes (avoid div-by-zero at BCs)
  n = length(u);
  ii = 2:n-1;
  nz = abs(ue(ii)) > 1e-15;
  if any(nz)
    max_rel = max(abs(u(ii(nz)) - ue(ii(nz))) ./ abs(ue(ii(nz))));
  else
    max_rel = max_abs;
  end
end


%% =====================================================================
%%  BDF3/EXT3 UNSTEADY SOLVER  (restricted system with boundary lifting)
%% =====================================================================
function u_final = bdf3ext3(Bbar, Cbar, Abar, n, dt, nsteps, nu, ...
                            u_left_fn, right_dirichlet, u_right_fn)
% Solves   Bbar u_t = -Cbar u - nu Abar u   (f = 0)
%
% IMEX splitting:
%   implicit — diffusion  (-nu Abar u)
%   explicit — advection  (-Cbar u)       extrapolated via EXT
%
% BDF time derivative on LHS, EXT extrapolation for advection on RHS.
% Boundary conditions are enforced with the restriction/prolongation
% operators described in adv_diff_hw4.pdf:
%   u = R' * u0 + ub,   H = R * Hbar * R'
%
% Startup: BDF1/EXT1 (step 1), BDF2/EXT2 (step 2), then BDF3/EXT3.

  % BDF coefficients  beta_0 / dt  for the LHS mass term
  b0 = [1, 3/2, 11/6];

  R = build_restriction(n, true, right_dirichlet);
  Rt = R';

  % Pre-build the three Helmholtz matrices and their restricted forms
  Hbar = cell(3, 1);
  H = cell(3, 1);
  for ord = 1:3
    Hbar{ord} = (b0(ord) / dt) * Bbar + nu * Abar;
    H{ord} = R * Hbar{ord} * Rt;
  end

  % Solution history:  uh{1} = u^{n},  uh{2} = u^{n-1},  uh{3} = u^{n-2}
  uh = { zeros(n, 1),  zeros(n, 1),  zeros(n, 1) };

  for step = 1:nsteps
    t_new = step * dt;

    if step == 1
      %  BDF1/EXT1
      bdf_rhs = uh{1};
      ext_rhs = uh{1};
      ord = 1;
    elseif step == 2
      %  BDF2/EXT2
      bdf_rhs = 2 * uh{1} - 0.5 * uh{2};
      ext_rhs = 2 * uh{1} -       uh{2};
      ord = 2;
    else
      %  BDF3/EXT3
      bdf_rhs = 3 * uh{1} - 1.5 * uh{2} + (1/3) * uh{3};
      ext_rhs = 3 * uh{1} - 3   * uh{2} +         uh{3};
      ord = 3;
    end

    ub = zeros(n, 1);
    ub(1) = u_left_fn(t_new);
    if right_dirichlet
      ub(n) = u_right_fn(t_new);
    end

    rhs = (1 / dt) * (Bbar * bdf_rhs) - Cbar * ext_rhs - Hbar{ord} * ub;

    u0_new = H{ord} \ (R * rhs);
    u_new = Rt * u0_new + ub;

    % Shift history
    uh{3} = uh{2};
    uh{2} = uh{1};
    uh{1} = u_new;
  end

  u_final = uh{1};
end


%% =====================================================================
%%  PLOTTING HELPERS
%% =====================================================================
function plot_q1(x, u, c, nu, f, L, max_err, ttl, fname)
  xf  = linspace(0, L, 10001)';
  uef = exact_steady(xf, c, nu, f, L);
  uhf = interp1(x, u, xf, 'linear');

  fig = figure('visible', 'off');

  subplot(2, 1, 1);
  plot(xf, uef, 'k-', 'linewidth', 2);  hold on;
  plot(x,  u,   'bo-', 'linewidth', 1.2, 'markersize', 4);  hold off;
  xlabel('x');  ylabel('u(x)');
  title(ttl);
  legend('Exact', 'FEM', 'location', 'northwest');
  grid on;

  subplot(2, 1, 2);
  plot(xf, abs(uhf - uef), 'r-', 'linewidth', 1.5);
  xlabel('x');  ylabel('|u_h - u_{exact}|');
  title(sprintf('Max abs. error = %.3e', max_err));
  grid on;

  print(fig, fname, '-dpng', '-r200');
  close(fig);
end


%% =====================================================================
%%  WRITE RESULTS FILE
%% =====================================================================
function write_results(abs_a, rel_a, abs_b, rel_b, s_b, ...
                       abs_c, rel_c, s_opt, ua, ub, x, c, nu, t_end)
  fid = fopen('hw3_results.txt', 'w');

  fprintf(fid, 'CS555 HW4 Results Summary\n');
  fprintf(fid, '=========================\n\n');

  fprintf(fid, 'Q1(a): Uniform mesh, E = 100\n');
  fprintf(fid, '  Max pointwise relative error: %.6e\n', rel_a);
  fprintf(fid, '  Max absolute error:           %.6e\n\n', abs_a);

  fprintf(fid, 'Q1(b): Geometric mesh, E = 20, s = %.2f\n', s_b);
  fprintf(fid, '  Max pointwise relative error: %.6e\n', rel_b);
  fprintf(fid, '  Max absolute error:           %.6e\n\n', abs_b);

  fprintf(fid, 'Q1(c): Optimal s for E = 20\n');
  fprintf(fid, '  Search:                       s in [0.50, 1.05] with ds = 0.005\n');
  fprintf(fid, '  Optimization metric:          max pointwise relative error\n');
  fprintf(fid, '  Best s:                       %.3f\n', s_opt);
  fprintf(fid, '  Max pointwise relative error: %.6e\n', rel_c);
  fprintf(fid, '  Max absolute error:           %.6e\n\n', abs_c);

  fprintf(fid, 'Q2(a): Dirichlet BCs, t = %.2f\n', t_end);
  fprintf(fid, '  min(u) = %+.6e,  max(u) = %+.6e,  u(1) = %+.6e\n\n', ...
          min(ua), max(ua), ua(end));

  fprintf(fid, 'Q2(b): Neumann at x=L, t = %.2f\n', t_end);
  fprintf(fid, '  min(u) = %+.6e,  max(u) = %+.6e,  u(1) = %+.6e\n', ...
          min(ub), max(ub), ub(end));

  fprintf(fid, '\n');
  write_problem2_discussion(fid, x, ua, ub, c, nu);

  fclose(fid);
end


%% =====================================================================
%%  RESTRICTION / PROLONGATION
%% =====================================================================
function R = build_restriction(n, left_dirichlet, right_dirichlet)
  active = true(1, n);
  if left_dirichlet
    active(1) = false;
  end
  if right_dirichlet
    active(n) = false;
  end

  free = find(active);
  R = sparse(1:length(free), free, 1, length(free), n);
end


%% =====================================================================
%%  WRITTEN DISCUSSION FOR PROBLEM 2
%% =====================================================================
function write_problem2_discussion(fid, x, ua, ub, c, nu)
  h = x(2) - x(1);
  pe_h = c * h / (2 * nu);
  tail_start = max(1, length(x) - 20);
  osc_a = count_sign_changes(diff(ua(tail_start:end)));
  osc_b = count_sign_changes(diff(ub(tail_start:end)));

  fprintf(fid, 'Discussion for Q2(a)\n');
  fprintf(fid, '  The Dirichlet outflow case shows strong nonphysical oscillations near x = 1.\n');
  fprintf(fid, '  In the last 20 intervals, the slope changes sign %d times, and the solution overshoots to %.6e before dropping to u(1,t) = 0.\n', ...
          osc_a, max(ua));
  fprintf(fid, '  This is consistent with the expected behavior of unstabilized linear FEM in an advection-dominated regime: h = %.4f gives cell Peclet number Pe_h = c h / (2 nu) = %.2f > 1.\n\n', ...
          h, pe_h);

  fprintf(fid, 'Discussion for Q2(b)\n');
  fprintf(fid, '  Replacing the right Dirichlet condition with the natural Neumann condition du/dx(1,t) = 0 removes the forced outflow value.\n');
  fprintf(fid, '  The final profile is much smoother and nearly monotone near the outflow, with %d slope reversals in the last 20 intervals and u(1,t) = %.6e.\n\n', ...
          osc_b, ub(end));

  fprintf(fid, 'Difference between Q2(a) and Q2(b)\n');
  fprintf(fid, '  Part (a) forces the advected solution to match u(1,t) = 0 at the outflow, which creates the large overshoot / undershoot pattern seen in the plot.\n');
  fprintf(fid, '  Part (b) behaves like a more natural outflow boundary, so the large oscillations are greatly reduced and the curve stays bounded near x = 1.\n');
end


%% =====================================================================
%%  SIGN-CHANGE COUNTER
%% =====================================================================
function n_changes = count_sign_changes(v)
  tol = 1e-12 * max(1, max(abs(v)));
  s = sign(v);
  s(abs(v) <= tol) = 0;
  s = s(s ~= 0);

  if length(s) < 2
    n_changes = 0;
  else
    n_changes = sum(s(2:end) ~= s(1:end-1));
  end
end
