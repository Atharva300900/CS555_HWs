% Plot script for 1c: BDF2 for both ICs

figure(1);

% --- BDF2 with u0 = sin(pi*x) ---
subplot(1,2,1);
N = 1500; nstep = 200;
n  = N-1; nu = .20; Lx = 1; T = 1.0;
dt = T/nstep; iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
H_eb  = I + dt*nu*A;
H_bdf = (3/2)*I + dt*nu*A;
sx = sin(pi*x);
u0 = sx; u = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',1.5); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep; time=k*dt;
  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',1.5); end
  if k==1; uold=u; u=H_eb\u;
  else rhs=2*u-(1/2)*uold; uold=u; u=H_bdf\rhs; end
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2);
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',14); ylabel('u(x,t)','FontSize',14);
title('BDF2,  u_0 = sin(\pi x)','fontsize',13);

% --- BDF2 with u0 = pi/4 ---
subplot(1,2,2);
N = 1500; nstep = 200;
n  = N-1; dt = T/nstep; iostep = floor(nstep/10);
L  = Lx; a_var;
I  = speye(n);
H_eb  = I + dt*nu*A;
H_bdf = (3/2)*I + dt*nu*A;
sx = sin(pi*x);
u0 = (pi/4)+0*x; u = u0;

ub=[0;u0;0]; hold off; plot(xb,ub,'k-','linewidth',1.5); hold on;
plot([0 1],[0 0],'k-');
for k=1:nstep; time=k*dt;
  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',1.5); end
  if k==1; uold=u; u=H_eb\u;
  else rhs=2*u-(1/2)*uold; uold=u; u=H_bdf\rhs; end
end
ub=[0;u;0]; plot(xb,ub,'b-','linewidth',2);
axis([0 1 -.5 1.3]); axis square;
xlabel('x','FontSize',14); ylabel('u(x,t)','FontSize',14);
title('BDF2,  u_0 = \pi/4','fontsize',13);

print('-dpng','-r150','plot_1c.png');
