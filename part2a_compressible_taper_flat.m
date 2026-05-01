% part2a_compressible_taper_flat.m
% CS555 HW6 Slider Bearing, Part 2a
% Compressible taper-flat slider using fixed-point iteration.

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

% Geometry-dependent height terms
he = profile_taper_flat(xe,L,T,Ta,h1,h2);
h3 = he.^3;

% Assemble Couette matrix C, same as in the incompressible code
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

% Compressible taper-flat case: p = 0 on the entire boundary
boundary_nodes = boundedges([x,y],t);
R = restriction(nb,boundary_nodes);
C = R*Cb*R';
rhs = R*(Cb*(patm*ones(nb,1)));

% Fixed-point iteration for pa = p + patm in the diffusion coefficient
Pb = zeros(nb,1);              % relative pressure, including boundary nodes
pa = patm*ones(nb,1);          % absolute pressure
max_iter = 50;
tol_iter = 1e-10;
rel_change = Inf;

for iter = 1:max_iter
    pae = sum(pa(t'),1)'/3;    % elementwise absolute pressure at centroids
    coef = (pae.*h3)/(12*mu);  % pa*h^3/(12 mu)
    coef = spdiags(coef,0,E,E);
    coef = kron(coef,speye(3));
    An = coef*AL;

    A = R*(Q'*An*Q)*R';
    P = (A - C) \ rhs;
    Pb_new = R'*P;

    rel_change = norm(Pb_new-Pb)/max(norm(Pb_new),1);
    Pb = Pb_new;
    pa = Pb + patm;

    if rel_change < tol_iter
        break;
    end
end

% Part 2a integrals. Use relative pressure p, not absolute pressure pa.
F  = ones(1,nb)*Bb*Pb;
xp = (x')*Bb*Pb/F;

fprintf('\nPart 2a: Compressible taper-flat case\n');
fprintf('Ex = %d, Ey = %d\n',Ex,Ey);
fprintf('h2 = %.8e m\n',h2);
fprintf('gamma = %.8e rad\n',gamma);
fprintf('Fixed-point iterations = %d\n',iter);
fprintf('Final relative pressure change = %.8e\n',rel_change);
fprintf('F = %.8e N\n',F);
fprintf('xp = %.8e m\n',xp);
fprintf('xp/L = %.8f\n',xp/L);
fprintf('Average relative pressure F/(L W) = %.8e Pa\n',F/(L*W));
fprintf('Average relative pressure F/(L W) = %.8f atm\n',F/(L*W)/patm);
fprintf('max relative pressure = %.8e Pa\n',max(Pb));
fprintf('max absolute pressure = %.8e Pa\n',max(pa));

% Pressure mesh plot using relative pressure p
Pmax = max(abs(Pb));
figure;
trimesh(t,x,y,W*Pb/Pmax);
axis([0 L -W/1.9 W/1.9 0 W]);
axis equal;
xlabel('x',fs,20);
ylabel('y',fs,20);
zlabel('W p / p_{max}',fs,20);
title('Part 2a: Compressible taper-flat pressure',fs,15);
