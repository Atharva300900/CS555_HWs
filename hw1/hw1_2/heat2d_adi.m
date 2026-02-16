
%
% HEAT EQUATION, u_t = nu \nabla^2 u, USING m-POINT FINITE DIFFERENCE
% CN with ADI (Alternating Direction Implicit)
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
%  Set up ADI tridiagonal operators (m x m each)
%
%  H_lx acts as left-multiply by T_lx:   T_lx * U
%  H_ly acts as right-multiply by T_ly:  U * T_ly
%  H_rx acts as left-multiply by T_rx:   T_rx * U
%  H_ry acts as right-multiply by T_ry:  U * T_ry
%

Tlx = Ix + (dt/2)*Ax;    % m x m tridiagonal
Trx = Ix - (dt/2)*Ax;
Tly = Iy + (dt/2)*Ay;    % m x m tridiagonal
Try_ = Iy - (dt/2)*Ay;   % Try is reserved word, use Try_

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

% Pre-factor the tridiagonal matrices (sparse \ is already efficient for tridiag)
% Octave/Matlab handles sparse tridiagonal solves in O(n)

for istep=1:nstep; time=istep*dt;

  % RHS: apply H_rx then H_ry
  %   Step R1: V = T_rx * U    (left multiply)
  %   Step R2: W = V * T_ry    (right multiply, T_ry symmetric)
  V = Trx * U;
  W = V * Try_;

  % Solve: H_lx then H_ly
  %   Step S1: T_lx * Z = W  -->  Z = T_lx \ W  (n_y tridiag solves of size n_x)
  Z = Tlx \ W;

  %   Step S2: Z * T_ly = U_new --> U_new' = T_ly \ Z'  (n_x tridiag solves of size n_y)
  U = (Tly \ Z')';

end;

Uex = exp(lamxy*time)*U0;
Err = Uex-U;
eo  = e2;
e2  = norm(Err,"fro") / norm(Uex,"fro");
ratio = eo/e2;

format shorte;
disp([N dt nstep time e2 ratio])
