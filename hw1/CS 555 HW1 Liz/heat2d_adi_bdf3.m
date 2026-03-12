function e2 = heat2d_adi_bdf3(N, dt, T)
% ADI-BDF3 for u_t = -K u, where K(U) = Ax*U + U*Ay'
% Bootstraps first two steps using ADI-CN.
% Uses linear extrapolation U_hat = 2U^{n-1} - U^{n-2} to correct ADI O(dt^2) term.

m = N - 1;
Lx = 1; Ly = 1;
nu = 0.05;

nstep = ceil(T/dt);
dt = T/nstep;

dx = Lx/N; xb = dx*(0:N)'; x = xb(2:end-1);
dy = Ly/N; yb = dy*(0:N)'; y = yb(2:end-1);

e = ones(m,1);
A = spdiags([-e 2*e -e], -1:1, m, m);

Ax = nu*A/(dx*dx);
Ay = nu*A/(dy*dy);
I  = speye(m);

[X,Y] = ndgrid(x,y);
kx = 1; ky = 3;
U0 = sin(kx*pi*X/Lx).*sin(ky*pi*Y/Ly);

thx  = pi*kx*dx/Lx;
thy  = pi*ky*dy/Ly;
lamx = 2*(1-cos(thx))/(dx*dx);
lamy = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx+lamy);

% BDF3 coefficients (divided by beta0)
beta0 = 11/6;
c1 = 18/11;
c2 = -9/11;
c3 = 2/11;
alpha = dt/beta0;

% ADI factors for BDF3 LHS in matrix form:
% (I + alpha*Ax) U (I + alpha*Ay)' = RHS + alpha^2*(Ax U Ay')
Tx = I + alpha*Ax;
Ty = I + alpha*Ay;

Lx_fac = chol(Tx,'lower');
Ly_fac = chol(Ty,'lower');

% Bootstrap with ADI-CN for first two steps:
TxL = I + (dt/2)*Ax;   TxR = I - (dt/2)*Ax;
TyL = I + (dt/2)*Ay;   TyR = I - (dt/2)*Ay;

Lx_cn = chol(TxL,'lower');
Ly_cn = chol(TyL,'lower');

% History: U^{0}, U^{1}, U^{2}
U_nm3 = U0;

% n=1 via ADI-CN
RHS = (TxR*U_nm3) * (TyR');
Wt  = Ly_cn'\(Ly_cn\RHS');
W   = Wt';
U_nm2 = Lx_cn'\(Lx_cn\W);

% n=2 via ADI-CN
RHS = (TxR*U_nm2) * (TyR');
Wt  = Ly_cn'\(Ly_cn\RHS');
W   = Wt';
U_nm1 = Lx_cn'\(Lx_cn\W);

time = 2*dt;

% BDF3 from n=3 onward
for n = 3:nstep
    time = n*dt;

    RHS0 = c1*U_nm1 + c2*U_nm2 + c3*U_nm3;

    % 2nd order surrogate for U^n
    U_hat = 2*U_nm1 - U_nm2;

    % Correct ADI splitting term: alpha^2 * Ax*U_hat*Ay'
    RHS = RHS0 + (alpha*alpha) * (Ax * U_hat * (Ay'));

    % Solve: Tx * U * Ty' = RHS
    % First solve U * Ty' = (Tx \ RHS)
    Z  = Lx_fac'\(Lx_fac\RHS);        % Z = Tx \ RHS
    Wt = Ly_fac'\(Ly_fac\Z');         % Wt = Ty \ Z'
    U  = Wt';

    % shift history
    U_nm3 = U_nm2;
    U_nm2 = U_nm1;
    U_nm1 = U;
end

Uex = exp(lamxy*time)*U0;
e2  = norm(Uex-U_nm1, "fro") / norm(Uex, "fro");
end