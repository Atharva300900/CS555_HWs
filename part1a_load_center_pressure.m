% part1a_load_center_pressure.m
% CS555 HW6 Slider Bearing, Part 1a
% Incompressible taper-flat slider: compute load F and center-of-pressure xp.

clear; close all; clc;
hdr;

Ex = 120;
Ey = 40;

rho = 1.225;         % kg/m^3      density of air
nu_air  = 1.5e-5;    % m^2/s       kinematic viscosity
mu  = rho*nu_air;    %             dynamic viscosity
U   = 20;            % m/s         speed of plate
L   = .0041;         % m           slider length
W   = .15625*L;      % m           slider width
T   = .09375*L;      % m           slider taper length
Ta  = pi/180;        % rad         taper angle, 1 degree for taper-flat case

h1  = 3.70e-7;       % m           leading edge gap, without taper
h2  = 2.50e-7;       % m           trailing edge gap
gamma = atan((h1-h2)/L); % rad     pitch angle
patm = 101325;       % N/m^2       atmospheric pressure

[x,y,t] = box_elem(Ex,Ey,L,W);
E  = 2*Ex*Ey;

nv = 3;              % linear triangles
xL = x(t');
yL = y(t');
xe = sum(xL,1)'/nv;  % centroids for coefficient quadrature

[AL,BL,Q,t,areaL] = abqfem([x y],t);
Bb = Q'*BL*Q;
nb = size(Q,2);

% Variable diffusion coefficient h^3/(12 mu)
he = profile_taper_flat(xe,L,T,Ta,h1,h2);
h3 = he.^3;
coef = (1/(12*mu))*h3;
coef = spdiags(coef,0,E,E);

Iv = speye(nv);
coef = kron(coef,Iv);
An = coef*AL;

% Assemble Couette forcing matrix C
order = 1;
nv = 3*order;

Nq = 3;
[z,w] = trigausspoints(Nq);
rq = z(:,1);
sq = z(:,2);
Bh = diag(w);

nq = length(rq);
Dr = zeros(nq,nv);
Ds = Dr;
Jq = zeros(nq,nv);
for k = 1:nv
    [Dr(:,k), Ds(:,k)] = basis_deriv_12(rq,sq,k,order);
    Jq(:,k) = basis_tri_12(rq,sq,k,order);
end

Xr  = Dr*xL;  Yr = Dr*yL;
Xs  = Ds*xL;  Ys = Ds*yL;
Jac = Xr.*Ys - Xs.*Yr;
Jmin = min(min(Jac));
if Jmin <= 0
    error('vanishing Jacobian');
end
Rx = Ys ./ Jac;
Sx = -Yr ./ Jac;
Bq = Bh*Jac;

hL = profile_taper_flat(xL,L,T,Ta,h1,h2);
hq = Jq*hL;
Uh = (U/2)*(Bq.*hq);

Ie = speye(E);
DLr = kron(Ie,Dr');
DLs = kron(Ie,Ds');
ULr = spdiags(reshape(Rx.*Uh,nq*E,1),0,nq*E,nq*E);
ULs = spdiags(reshape(Sx.*Uh,nq*E,1),0,nq*E,nq*E);
JL  = kron(Ie,Jq);

CL = (DLr*ULr + DLs*ULs)*JL;
Cb = Q'*CL*Q;

% Original taper-flat case: Dirichlet pressure on the entire boundary
boundary_nodes = boundedges([x,y],t);
R = restriction(nb,boundary_nodes);

A = R*(Q'*An*Q)*R';
rhs = R*(Cb*ones(nb,1));

P  = A \ rhs;
Pb = R'*P;

% Part 1a integrals
F  = ones(1,nb)*Bb*Pb;
xp = (x')*Bb*Pb/F;

fprintf('\nPart 1a: Incompressible taper-flat case\n');
fprintf('Ex = %d, Ey = %d\n',Ex,Ey);
fprintf('h2 = %.8e m\n',h2);
fprintf('gamma = %.8e rad\n',gamma);
fprintf('F = %.8e N\n',F);
fprintf('xp = %.8e m\n',xp);
fprintf('xp/L = %.8f\n',xp/L);
fprintf('Average pressure F/(L W) = %.8e Pa\n',F/(L*W));
fprintf('Average pressure F/(L W) = %.8f atm\n',F/(L*W)/patm);

% Pressure mesh plot
Pmax = max(abs(Pb));
figure;
trimesh(t,x,y,W*Pb/Pmax);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
xlabel('x',fs,20);
ylabel('y',fs,20);
zlabel('W p / p_{max}',fs,20);
title('Part 1a: Incompressible taper-flat pressure',fs,15);
