function [solver_cfl, target_dt] = hw2_effective_cfl(N, grid_type, target_dt)
% Map the homework timestep choice onto burgers_solve's CFL input.
if nargin < 3 || isempty(target_dt)
    target_dt = 0.1 / N;
end
solver_cfl = target_dt * N;
end
