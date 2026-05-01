% part2bc_compare_cases.m
% CS555 HW6 Slider Bearing, Parts 2b and 2c
% Compare pressure distributions, load, and center-of-pressure for:
%   Case 1a: incompressible taper-flat
%   Case 1b: incompressible wedge verification configuration
%   Case 2a: compressible taper-flat

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
h1  = 3.70e-7;       % m           leading edge gap, without taper
h2  = 2.50e-7;       % m           trailing edge gap
gamma = atan((h1-h2)/L); % rad     pitch angle
patm = 101325;       % N/m^2       atmospheric pressure

% -------------------------
% Case definitions
% -------------------------
case_name = {'1a: incompressible taper-flat', ...
             '1b: incompressible wedge', ...
             '2a: compressible taper-flat'};

Ta_list = [pi/180, 0, pi/180];
use_compressible = [false, false, true];
use_wedge_bc = [false, true, false];

F = zeros(3,1);
xp = zeros(3,1);
Pmax = zeros(3,1);
Pmin = zeros(3,1);
iters = zeros(3,1);

Pstore = cell(3,1);
xstore = cell(3,1);
ystore = cell(3,1);
tstore = cell(3,1);

for c = 1:3
    Ta = Ta_list(c);

    [x,y,t] = box_elem(Ex,Ey,L,W);
    if ishandle(gcf)
        close(gcf); % box_elem makes a mesh plot
    end
    E  = 2*Ex*Ey;

    nv = 3;
    xL = x(t');
    yL = y(t');
    xe = sum(xL,1)'/nv;

    [AL,BL,Q,t,areaL] = abqfem([x y],t);
    Bb = Q'*BL*Q;
    nb = size(Q,2);

    he = profile_taper_flat(xe,L,T,Ta,h1,h2);
    h3 = he.^3;

    % Assemble Couette matrix C
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

    if use_wedge_bc(c)
        % Wedge verification: Dirichlet only at x = 0 and x = L.
        % Top and bottom are natural Neumann boundaries.
        tol_bc = 1e-12;
        boundary_nodes = find(abs(x) < tol_bc | abs(x-L) < tol_bc);
    else
        % Taper-flat cases: p = 0 on the full boundary.
        boundary_nodes = boundedges([x,y],t);
    end
    R = restriction(nb,boundary_nodes);

    if ~use_compressible(c)
        % Incompressible system: A p = C 1
        coef = (1/(12*mu))*h3;
        coef = spdiags(coef,0,E,E);
        coef = kron(coef,speye(3));
        An = coef*AL;

        A = R*(Q'*An*Q)*R';
        rhs = R*(Cb*ones(nb,1));
        P = A \ rhs;
        Pb = R'*P;
        iters(c) = 1;
    else
        % Compressible fixed-point system: (A(pa) - C) p = C patm
        C = R*Cb*R';
        rhs = R*(Cb*(patm*ones(nb,1)));

        Pb = zeros(nb,1);
        pa = patm*ones(nb,1);
        max_iter = 50;
        tol_iter = 1e-10;

        for iter = 1:max_iter
            pae = sum(pa(t'),1)'/3;
            coef = (pae.*h3)/(12*mu);
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
        iters(c) = iter;
    end

    % Load and center-of-pressure use relative pressure p for all cases.
    F(c)  = ones(1,nb)*Bb*Pb;
    xp(c) = (x')*Bb*Pb/F(c);
    Pmax(c) = max(Pb);
    Pmin(c) = min(Pb);

    Pstore{c} = Pb;
    xstore{c} = x;
    ystore{c} = y;
    tstore{c} = t;
end

% -------------------------
% Part 2b: vertically stacked pressure distributions
% -------------------------
figure;

for c = 1:3
    subplot(3,1,c);

    Pb = Pstore{c};
    x  = xstore{c};
    y  = ystore{c};
    t  = tstore{c};

    scale = max(abs(Pb));

    trimesh(t,x,y,W*Pb/scale);

    axis([0 L -W/1.9 W/1.9 -W W]);
    axis equal;

    xlabel('x');
    ylabel('y');
    zlabel('W p / max|p|');
    title(case_name{c},'fontsize',11);

    view(35,25);
end

sgtitle('Pressure distributions for Parts 1a, 1b, and 2a');

% -------------------------
% Part 2c: table
% -------------------------
fprintf('\nParts 2b and 2c: comparison table\n');
fprintf('Ex = %d, Ey = %d\n',Ex,Ey);
fprintf('h2 = %.8e m, gamma = %.8e rad\n\n',h2,gamma);
fprintf('%34s %14s %14s %14s %14s %8s\n', ...
    'case','F [N]','xp/L','p_min [Pa]','p_max [Pa]','iters');
for c = 1:3
    fprintf('%34s %14.8e %14.8f %14.8e %14.8e %8d\n', ...
        case_name{c},F(c),xp(c)/L,Pmin(c),Pmax(c),iters(c));
end

fprintf('\nNote: F and xp are computed using the relative pressure p.\n');
fprintf('For the compressible case, pa = p + patm, so integrating pa would add an artificial ambient-pressure load.\n');
