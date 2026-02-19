
%
% HEAT EQUATION, u_t = nu \nabla^2 u, USING m-POINT FINITE DIFFERENCE
% CN with permuted Cholesky factorization (symamd)
%
% For CN we must also multiply by HR each step. We permute HR as well
% and keep u in permuted order throughout, permuting/unpermuting only
% at the start and end (equations 18-24 adapted for CN).
%

hdr

if ~exist('N','var'); N=150; end;
m=N-1;

Lx=1; Ly=1;
T = 1.2; nstep = ceil(T/dt); dt=T/nstep;
nu = .05;

dx=Lx/N; xb = dx*[0:N]'; x=xb(2:end-1);
e = ones(m,1); A = spdiags([-e 2*e -e],-1:1, m,m);
Ax = nu*A/(dx*dx);
Ix = speye(m);


dy=Ly/N; yb = dy*[0:N]'; y=yb(2:end-1);
e = ones(m,1); A = spdiags([-e 2*e -e],-1:1, m,m);
Ay = nu*A/(dy*dy);
Iy = speye(m);

[X,Y]=ndgrid(x,y);
[Xb,Yb]=ndgrid(xb,yb);

%
%  Set up standard CN operators
%

HL = kron(Iy,Ix) + (dt/2)*(kron(Iy,Ax)+kron(Ay,Ix));
HR = kron(Iy,Ix) - (dt/2)*(kron(Iy,Ax)+kron(Ay,Ix));

%
%  Set up Exact Solution + RHS
%

kx = 1; ky = 3;

U0 = sin(kx*pi*X/Lx).*sin(ky*pi*Y/Ly);
U  = U0;

lamxy = -nu*( (kx*pi/Lx)^2 + (ky*pi/Ly)^2 );
thx   = pi*kx*dx/Lx;
thy   = pi*ky*dy/Ly;
lamx  = 2*(1-cos(thx))/(dx*dx);
lamy  = 2*(1-cos(thy))/(dy*dy);
lamxy = -nu*(lamx+lamy);              % Judge accuracy by discrete lambdas

u = reshape(U,m*m,1);

% Permuted Cholesky factorization
p    = symamd(HL);
L    = chol(HL(p,p),'lower');
%HRp  = HR(p,p);                  % pre-permute HR once outside the loop

%u = u(p);                        % permute u once before the loop


for istep=1:nstep; time=istep*dt;
  rhs = HR*u;
  rhs_p = rhs(p);
  u_p = L'\( L \ (rhs_p) );       % everything stays in permuted space
  u = zeros(size(u));
  u(p) = u_p;
end
                       % unpermute u once after the loop
U = reshape(u,m,m);

Uex = exp(lamxy*time)*U0;
Err = Uex-U;
e2  = norm(Err,"fro") / norm(Uex,"fro");
eo  = e2;
ratio = eo/e2;

format shorte;
disp([N dt nstep time e2 ratio])
