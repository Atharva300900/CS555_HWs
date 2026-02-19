figure(1);

% --- EB plot ---
subplot(1,2,1);
N = 1500; nstep = 200;
n  = N-1; nu = .20; Lx = 1; T = 1.0;
dt = T/nstep; iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
H  = I + dt*nu*A;
sx = sin(pi*x);
u0 = sx; u = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',2.5); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep; time=k*dt;
  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2.5); end
  u = H\u;
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2.5);
ub=[0;u0;0]; plot(xb,ub,'k-','linewidth',2.5);
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',20); ylabel('u(x,t)','FontSize',20);
title('EB,  u_0 = sin(\pi x)','fontsize',30);

% --- CN plot ---
subplot(1,2,2);
N = 1500; nstep = 200;
n  = N-1; dt = T/nstep; iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
HL = I + (dt/2)*nu*A;
HR = I - (dt/2)*nu*A;
sx = sin(pi*x);
u0 = sx; u = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',2.5); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep; time=k*dt;
  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2.5); end
  u = HL\(HR*u);
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2.5);
ub=[0;u0;0]; plot(xb,ub,'k-','linewidth',2.5);
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',20); ylabel('u(x,t)','FontSize',20);
title('CN,  u_0 = sin(\pi x)','fontsize',30);

print('-dpng','-r150','plot_1a.png');
