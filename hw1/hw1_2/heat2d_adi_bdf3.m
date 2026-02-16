
%
% HEAT EQUATION, u_t = nu \nabla^2 u, USING m-POINT FINITE DIFFERENCE
% ADI-BDF3 (bootstrapped with 2 steps of ADI-CN)
%
% The ADI splitting for BDF3 introduces an O(dt^3) error per step:
%   H_ADI = H_BDF3 + (6*dt^2/11)*(Ay kron Ax)
% To maintain 3rd-order accuracy, we correct the RHS using a linear
% extrapolation u_ext = 2*u^n - u^{n-1} (2nd-order accurate):
%   rhs_corrected = rhs_bdf3 + (6*dt^2/11)*(Ax * U_ext * Ay)
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
%  Set up ADI operators
%

% CN operators (for bootstrap)
Tlx_cn  = Ix + (dt/2)*Ax;
Trx_cn  = Ix - (dt/2)*Ax;
Tly_cn  = Iy + (dt/2)*Ay;
Try_cn  = Iy - (dt/2)*Ay;

% BDF3 operators:  H_BDF3 = (11/6)I + dt*A
% ADI factored:    (11/6)(I + c*Ax_full)(I + c*Ay_full), c = 6*dt/11
c_bdf = 6*dt/11;
Tlx_bdf = Ix + c_bdf*Ax;
Tly_bdf = Iy + c_bdf*Ay;

% Splitting correction coefficient
split_coeff = 6*dt^2/11;

%
%  Set up Exact Solution
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

%
%  Time stepping
%

% Store history as m x m matrices
Unm2 = U;   % u^{n-2}
Unm1 = U;   % u^{n-1}
Un   = U;   % u^{n}

for istep=1:nstep; time=istep*dt;

  if istep <= 2
    % Bootstrap: ADI-CN
    %   RHS: V = Trx * U, W = V * Try
    %   Solve: Z = Tlx \ W,  U_new = (Tly \ Z')'
    V = Trx_cn * Un;
    W = V * Try_cn;
    Z = Tlx_cn \ W;
    Unew = (Tly_cn \ Z')';

  else
    % BDF3: (11/6)u^{n+1} + dt*A*u^{n+1} = 3u^n - (3/2)u^{n-1} + (1/3)u^{n-2}
    %
    % ADI factored solve with splitting correction:
    %   rhs = 3*Un - (3/2)*Unm1 + (1/3)*Unm2 + split_coeff*(Ax*Uext*Ay)
    %   where Uext = 2*Un - Unm1  (linear extrapolation, O(dt^2) accurate)
    %
    %   (11/6) * Tlx_bdf * U^{n+1} * Tly_bdf = rhs
    %   Step 1: Tlx_bdf * W = (6/11) * rhs
    %   Step 2: U^{n+1} = (Tly_bdf \ W')'

    Uext = 2*Un - Unm1;
    RHS  = 3*Un - (3/2)*Unm1 + (1/3)*Unm2 + split_coeff*(Ax*Uext*Ay);

    W    = Tlx_bdf \ ((6/11)*RHS);
    Unew = (Tly_bdf \ W')';

  end

  % Shift history
  Unm2 = Unm1;
  Unm1 = Un;
  Un   = Unew;

end;

U = Un;
Uex = exp(lamxy*time)*U0;
Err = Uex-U;
eo  = e2;
e2  = norm(Err,"fro") / norm(Uex,"fro");
ratio = eo/e2;

format shorte;
disp([N dt nstep time e2 ratio])
