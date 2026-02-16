
% Extra Credit: Estimate the magnitude of the next leading-order term
% in the decaying square wave at t=1.
%
% u0 = pi/4 has Fourier sine series:
%   u0 = sum_{k=1,3,5,...} b_k sin(k*pi*x),   b_k = 1/k
%
% Exact solution:
%   u(x,t) = sum_{k=1,3,5,...} (1/k) exp(-nu*(k*pi)^2 * t) sin(k*pi*x)
%
% At t=1, the solution is approximately the k=1 eigenmode.
% The next leading-order term is k=3.

nu = 0.2;
t  = 1.0;

fprintf('\n===== Extra Credit: Square wave Fourier mode magnitudes at t=1 =====\n\n');
fprintf('  k     b_k      exp(-nu*(k*pi)^2*t)     |coefficient|\n');
fprintf('  ---   ------   ---------------------   ----------------\n');

for k = 1:2:11
  bk   = 1/k;
  decay = exp(-nu*(k*pi)^2*t);
  coeff = bk * decay;
  fprintf('  %2d    1/%d      %14.6e          %14.6e\n', k, k, decay, coeff);
end

% The k=1 and k=3 terms specifically:
c1 = exp(-nu*(pi)^2*t);
c3 = (1/3)*exp(-nu*(3*pi)^2*t);

fprintf('\n--- Key comparison ---\n');
fprintf('k=1 coefficient:  %e\n', c1);
fprintf('k=3 coefficient:  %e  (next leading-order term)\n', c3);
fprintf('Ratio |c3|/|c1|:  %e\n\n', abs(c3/c1));

% Compare to FD errors from preceding questions at various nstep values
fprintf('--- Comparison to FD temporal errors at t=1 ---\n');
fprintf('  (from nstep=1600, N=1500 runs)\n\n');
fprintf('  EB  error (sin IC, nstep=1600):   ~1.22e-03   (1st order)\n');
fprintf('  CN  error (pi/4 IC, nstep=1600):  ~8.49e-01   (oscillations!)\n');
fprintf('  BDF2 error (sin IC, nstep=1600):  ~2.85e-07   (2nd order)\n');
fprintf('  BDF2 error (pi/4 IC, nstep=1600): ~2.41e-07   (2nd order)\n');
fprintf('\n');
fprintf('  k=3 term magnitude:               ~%8.2e\n\n', c3);

fprintf('Conclusion:\n');
fprintf('  The k=3 term at t=1 has magnitude ~%.2e,\n', c3);
fprintf('  which is FAR smaller than the EB and BDF2 temporal errors\n');
fprintf('  at all dt values tested. This confirms that the assumption\n');
fprintf('  that u(x,1) is approximately the k=1 eigenmode is excellent.\n');
fprintf('  The higher Fourier modes have decayed to negligible levels,\n');
fprintf('  and the FD temporal discretization error completely dominates\n');
fprintf('  over the Fourier truncation error.\n\n');

% Also compute L2 norm of the k=3 contribution relative to k=1
% L2 norm of sin(k*pi*x) on [0,1] = 1/sqrt(2)
% So L2 of k=3 term / L2 of k=1 term = |c3/c1|
fprintf('  Relative L2 contribution of k=3 vs k=1: %.2e\n\n', abs(c3/c1));
