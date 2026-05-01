% part1b_wedge_convergence.m
% CS555 HW6 Slider Bearing, Part 1b
% Wedge-bearing verification against the analytical 1D solution.

clear; close all; clc;
hdr;

Ex_list = [10 20 40 80 120 160];
Ey = 40;

rho = 1.225;         % kg/m^3      density of air
nu_air  = 1.5e-5;    % m^2/s       kinematic viscosity
mu  = rho*nu_air;    %             dynamic viscosity
U   = 20;            % m/s         speed of plate
L   = .0041;         % m           slider length
W   = .15625*L;      % m           slider width
T   = .09375*L;      % m           slider taper length
Ta  = 0;             % rad         no leading taper for wedge verification

h1  = 3.70e-7;       % m           leading edge gap, without taper
h2  = 2.50e-7;       % m           trailing edge gap

err = zeros(size(Ex_list));
F_num = zeros(size(Ex_list));
xp_num = zeros(size(Ex_list));

alpha  = h1/h2;
Lambda = 6*mu*U*L/(h2^2);
F_exact = Lambda*L*W/(1-alpha)^2*(log(alpha) + 2*(1-alpha)/(1+alpha));

for m = 1:length(Ex_list)
    Ex = Ex_list(m);

    [x,y,t] = box_elem(Ex,Ey,L,W);
    close(gcf); % box_elem makes a mesh plot, close it during the convergence loop
    E = 2*Ex*Ey;

    nv = 3;
    xL = x(t');
    yL = y(t');
    xe = sum(xL,1)'/nv;

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

    % Wedge verification boundary conditions:
    % Dirichlet only at x = 0 and x = L.
    % Top and bottom y boundaries are left free, giving natural Neumann conditions.
    tol = 1e-12;
    boundary_nodes = find(abs(x) < tol | abs(x-L) < tol);
    R = restriction(nb,boundary_nodes);

    A = R*(Q'*An*Q)*R';
    rhs = R*(Cb*ones(nb,1));

    P  = A \ rhs;
    Pb = R'*P;

    % Analytical wedge pressure at every mesh node
    H = alpha + (1-alpha)*(x/L);
    P_exact = alpha*Lambda/(1-alpha^2)*(1./H.^2 - 1/alpha^2) ...
            - Lambda/(1-alpha)*(1./H - 1/alpha);

    err(m) = max(abs(Pb - P_exact));
    F_num(m) = ones(1,nb)*Bb*Pb;
    xp_num(m) = (x')*Bb*Pb/F_num(m);
end

% Convergence plot
figure;
loglog(Ex_list,err,'o-',lw,2); hold on;
loglog(Ex_list,err(1)*(Ex_list(1)./Ex_list).^2,'k--',lw,1.5);
grid on;
xlabel('E_x',fs,16);
ylabel('max |p_h - p_{exact}|',fs,16);
title('Part 1b: Wedge-bearing convergence',fs,15);
legend('numerical error','O(E_x^{-2}) reference','location','southwest');

fprintf('\nPart 1b: Wedge-bearing convergence check\n');
fprintf('Ey = %d, Ta = 0, Neumann on top and bottom\n',Ey);
fprintf('Analytical wedge load F_exact = %.8e N\n\n',F_exact);
fprintf('%8s %18s %18s %18s %18s\n','Ex','max error','F_num','rel F error','xp/L');
for m = 1:length(Ex_list)
    rel_F_error = abs(F_num(m)-F_exact)/abs(F_exact);
    fprintf('%8d %18.8e %18.8e %18.8e %18.8f\n', ...
        Ex_list(m),err(m),F_num(m),rel_F_error,xp_num(m)/L);
end

% Optional: plot the finest-grid numerical and analytical pressure along y = 0
Ex = Ex_list(end);
% The arrays from the last loop correspond to the finest grid.
mid_nodes = find(abs(y) == min(abs(y)));
[~,ord] = sort(x(mid_nodes));
mid_nodes = mid_nodes(ord);
figure;
plot(x(mid_nodes),Pb(mid_nodes),'o-',lw,1.5); hold on;
plot(x(mid_nodes),P_exact(mid_nodes),'k--',lw,1.5);
grid on;
xlabel('x',fs,16);
ylabel('p(x)',fs,16);
title('Finest-grid wedge pressure comparison',fs,15);
legend('FEM','analytical','location','best');
