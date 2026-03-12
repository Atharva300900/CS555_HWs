function results = burgers_solve(N, nu, T, grid_type, form, CFL, t_save)
% BURGERS_SOLVE  Solve 1D Burgers equation on [0,1] using BDF3/EXT3

if nargin < 7
    t_save = [];
end

% Build operators
[A, D, M, x_int, xb] = build_operators(N, nu, grid_type);
n = length(x_int);

% Initial condition
u = sin(pi * x_int);

% The homework uses dt = CFL / N with |c| = 1 and dx = 1 / N.
dt = CFL / N;

% Adjust dt so final time is exactly T
nsteps = ceil(T / dt);
dt = T / nsteps;
t_all = (0:nsteps) * dt;

% Storage for s(t)
s_all = zeros(1, nsteps + 1);

% Initial slope measure
du = D * u;
s_all(1) = max(abs(du));

% Snapshot storage (store full solution including boundaries)
if ~isempty(t_save)
    save_steps = round(t_save / dt);
    if any(abs(save_steps*dt - t_save) > 1e-10)
        warning('Some t_save values are not exact multiples of dt; using nearest time step.');
    end
    u_snap = zeros(N+1, length(t_save));
    t_snap = save_steps * dt;

    for si = 1:length(t_save)
        if save_steps(si) == 0
            u_snap(:, si) = [0; u; 0];
        end
    end
else
    u_snap = [];
    t_snap = [];
    save_steps = [];
end

% History values
u_nm1 = zeros(n, 1);
u_nm2 = zeros(n, 1);

% Initial nonlinear term
if strcmpi(form, 'convective')
    f_n = u .* (D * u);
elseif strcmpi(form, 'conservation')
    f_n = 0.5 * (D * (u.^2));
else
    error('form must be ''convective'' or ''conservation''.');
end

f_nm1 = zeros(n, 1);
f_nm2 = zeros(n, 1);

% Precompute Helmholtz operators
H1 = 1.0   * M + dt * A;
H2 = 1.5   * M + dt * A;
H3 = (11/6)* M + dt * A;

for step = 1:nsteps

    if step == 1
        b0=1;    b1=-1;   b2=0;    b3=0;
        a1=1;    a2=0;    a3=0;
        H = H1;
    elseif step == 2
        b0=1.5;  b1=-2;   b2=0.5;  b3=0;
        a1=2;    a2=-1;   a3=0;
        H = H2;
    else
        b0=11/6; b1=-3;   b2=1.5;  b3=-1/3;
        a1=3;    a2=-3;   a3=1;
        H = H3;
    end

    rhs = -(b1*(M*u) + b2*(M*u_nm1) + b3*(M*u_nm2)) ...
          - dt*(a1*f_n + a2*f_nm1 + a3*f_nm2);

    u_new = H \ rhs;

    % Shift history
    u_nm2 = u_nm1;
    u_nm1 = u;
    f_nm2 = f_nm1;
    f_nm1 = f_n;

    u = u_new;

    % Recompute nonlinear term
    if strcmpi(form, 'convective')
        f_n = u .* (D * u);
    else
        f_n = 0.5 * (D * (u.^2));
    end

    % Diagnostic
    du = D * u;
    s_all(step + 1) = max(abs(du));

    % Save snapshots
    if ~isempty(save_steps)
        idx = find(save_steps == step);
        if ~isempty(idx)
            for k = 1:length(idx)
                u_snap(:, idx(k)) = [0; u; 0];
            end
        end
    end
end

[s_star, idx_star] = max(s_all);
t_star = t_all(idx_star);

results.s_star = s_star;
results.t_star = t_star;
results.nsteps = nsteps;
results.s_all  = s_all;
results.t_all  = t_all;
results.x_int  = x_int;
results.xb     = xb;
results.dt     = dt;
results.u_snap = u_snap;
results.t_snap = t_snap;
end
