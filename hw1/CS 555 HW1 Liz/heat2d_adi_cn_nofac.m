function e2 = heat2d_adi_cn_nofac(N, dt, T)
% ADI-CN without explicit factorization:
% (I + dt/2*Ay)(I + dt/2*Ax) U^n = (I - dt/2*Ay)(I - dt/2*Ax) U^{n-1}
% dropping O(dt^3) term

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
U  = U0;

thx  = pi*kx*dx/Lx;
thy  = pi*ky*dy/Ly;
lamx = 2*(1-cos(thx))/(dx*dx);
lamy = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx+lamy);

TxL = I + (dt/2)*Ax;   TxR = I - (dt/2)*Ax;
TyL = I + (dt/2)*Ay;   TyR = I - (dt/2)*Ay;

time = 0;
for istep = 1:nstep
    time = istep*dt;

    RHS = (TxR*U) * (TyR');

    % Solve (TyL) on the right via transpose trick:
    % Want W * TyL' = RHS  ->  (TyL * W') = RHS'
    Wt = TyL \ RHS';
    W  = Wt';

    % Solve (TxL) on the left:
    U = TxL \ W;
end

Uex = exp(lamxy*time)*U0;
e2  = norm(Uex-U, "fro") / norm(Uex, "fro");
end