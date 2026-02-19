function [e2, t_fac_out] = heat2d_adi_2(N, dt, T)
% heat2d_adi.m
% ADI version of CN for 2D heat equation.
% Uses the splitting:
%   (TyL ⊗ Ix)(Iy ⊗ TxL) u^n = (TyR ⊗ Ix)(Iy ⊗ TxR) u^{n-1}
% which corresponds to dropping the O(dt^3) coupling term.
hdr
if ~exist('N','var'); N=150; end;
m = N - 1;
Lx = 1; Ly = 1;
nu = 0.05;
Lx=1; Ly=1;
dt = 0.01;
T  = 1.2;

nu = .05;
nstep = ceil(T / dt);

dx = Lx / N; xb = dx * (0:N)'; x = xb(2:end-1);
dy = Ly / N; yb = dy * (0:N)'; y = yb(2:end-1);

e = ones(m,1);
A = spdiags([-e 2*e -e], -1:1, m, m);

Ax = nu * A / (dx*dx);
Ay = nu * A / (dy*dy);

Ix = speye(m);

[X, Y] = ndgrid(x, y);

% Build 1D Helmholtz-like operators
TxL = Ix + (dt/2) * Ax;
TxR = Ix - (dt/2) * Ax;

TyL = Ix + (dt/2) * Ay;
TyR = Ix - (dt/2) * Ay;

% Exact eigenmode IC
kx = 1; ky = 3;
U0 = sin(kx*pi*X/Lx) .* sin(ky*pi*Y/Ly);
U  = U0;

thx  = pi*kx*dx/Lx;
thy  = pi*ky*dy/Ly;
lamx = 2*(1-cos(thx))/(dx*dx);
lamy = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx + lamy);

% Factorizations (both SPD)
t_fac = tic;
Lx_fac = chol(TxL, 'lower');
Ly_fac = chol(TyL, 'lower');
t_fac_out = toc(t_fac);

time = 0;
for istep = 1:nstep
    time = istep * dt;

    % RHS = Hry * Hrx * u^{n-1}
    % In matrix form with U (m x m):
    %   Hrx: vec(TxR * U)
    %   Hry: vec((TxR * U) * TyR')
    RHS = (TxR * U) * (TyR');  % TyR is symmetric, but keep transpose for clarity

    % Solve (TyL ⊗ Ix) W = RHS  ->  W * TyL' = RHS  ->  W = (TyL \ RHS')'
    Wt = Ly_fac' \ (Ly_fac \ RHS'); % solves TyL * Wt = RHS'
    W  = Wt';

    % Solve (Iy ⊗ TxL) U = W  ->  TxL * U = W
    U = Lx_fac' \ (Lx_fac \ W);
end

Uex = exp(lamxy * time) * U0;

Err = Uex - U;
e2 = norm(Err, "fro") / norm(Uex, "fro");
format shorte;
disp([N dt nstep time e2])
end


