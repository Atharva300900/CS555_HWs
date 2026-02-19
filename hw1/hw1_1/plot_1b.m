% Plot script for 1b: u0 = pi/4, EB vs CN side by side

figure(1);

% --- EB plot ---
subplot(1,2,1);
N = 1500; nstep = 200;
n  = N-1;
nu = .20;
Lx = 1;
T  = 1.0;
dt = T/nstep;
iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
H  = I + dt*nu*A;
sx = sin(pi*x);
u0 = (pi/4)+0*x;
u  = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',2.5); hold on;
%ub=[0;u0;0]; hold off; fill([xb; flipud(xb)],[ub; zeros(size(ub))],'k'); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep;  time=k*dt;
  if mod(k-1,iostep)==0
    ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2.5);
  end
  u = H\u;
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2.5);
ub=[0;u0;0]; plot(xb,ub,'k-','linewidth',2.5);   % re-plot IC on top
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',20); ylabel('u(x,t)','FontSize',20);
title('EB,  u_0 = \pi/4,  N=1500,  nstep=200','fontsize',30);

% --- CN plot ---
subplot(1,2,2);
N = 1500; nstep = 200;
n  = N-1;
dt = T/nstep;
iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
HL = I + (dt/2)*nu*A;
HR = I - (dt/2)*nu*A;
sx = sin(pi*x);
u0 = (pi/4)+0*x;
u  = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',2.5); hold on;
%ub=[0;u0;0]; hold off; fill([xb; flipud(xb)],[ub; zeros(size(ub))],'k'); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep;  time=k*dt;
  if mod(k-1,iostep)==0
    ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2.5);
  end
  u = HL\(HR*u);
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2.5);
ub=[0;u0;0]; plot(xb,ub,'k-','linewidth',2.5);   % re-plot IC on top
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',20); ylabel('u(x,t)','FontSize',20);
title('CN,  u_0 = \pi/4,  N=1500,  nstep=200','fontsize',30);

set(gcf, 'Position', [100 100 1200 500]);
set(gcf, 'PaperPositionMode', 'auto');
