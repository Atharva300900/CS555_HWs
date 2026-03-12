function e2 = heat2d_cn(N, dt, T)
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

Ix = speye(m); Iy = speye(m);
[X,Y] = ndgrid(x,y);

K  = kron(Iy,Ax) + kron(Ay,Ix);
HL = kron(Iy,Ix) + (dt/2)*K;
HR = kron(Iy,Ix) - (dt/2)*K;

kx = 1; ky = 3;
U0 = sin(kx*pi*X/Lx).*sin(ky*pi*Y/Ly);
u  = reshape(U0, m*m, 1);

thx  = pi*kx*dx/Lx;
thy  = pi*ky*dy/Ly;
lamx = 2*(1-cos(thx))/(dx*dx);
lamy = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx+lamy);

L = chol(HL,'lower');

time = 0;
for istep = 1:nstep
    time = istep*dt;
    u = L'\(L\(HR*u));
end

U   = reshape(u, m, m);
Uex = exp(lamxy*time)*U0;
e2  = norm(Uex-U, "fro") / norm(Uex, "fro");
end