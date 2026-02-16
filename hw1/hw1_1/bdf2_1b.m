
%
% HEAT EQUATION, u_t = nu u_xx, USING n-POINT FINITE DIFFERENCE
% BDF2 time integration (bootstrapped with EB for first step)
% IC: u0 = pi/4
%

n  = N-1;          % Number of points
nu = .20;
Lx = 1;

T  = 1.0;          % Final time
dt = T/nstep;
iostep = floor(nstep/10);

L  = Lx; a_var;  %% Set A based on variable (here, Chebyshev) spacing

I  = speye(n);
H_eb  = I + dt*nu*A;            % EB matrix for bootstrap step
H_bdf = (3/2)*I + dt*nu*A;     % BDF2 left-hand side

t1 = pi/N;
sx = sin(pi*x);   lam_ex = -nu*pi*pi;

u0 = sx;          % Smooth IC
u0 = (pi/4)+0*x;  % Constant IC
u  = u0;

% PLOT AXIS and Initial Condition:
ub=[0;u0;0]; hold off; plot(xb,ub,'k-',[0 1],[0 0],'k-'); hold on;
axis([0 1 -.5 1.3]); axis square;
xlabel('-- x --','FontSize',14); ylabel('u(x,t)','FontSize',14);

for k=1:nstep;  time=k*dt;

  if mod(k-1,iostep)==0; ub=[0;u;0]; plot(xb,ub,'r-','linewidth',2);
   title(['BDF2, N = ' int2str(N) ',  Step = ' int2str([k])],'fontsize',15); pause(.1); end;

  if k==1
    % Bootstrap: Euler Backward for first step
    uold = u;
    u    = H_eb\u;                 name='BDF2';
  else
    % BDF2: (3/2)u^{l+1} = 2*u^l - (1/2)*u^{l-1} + dt*nu*A*u^{l+1}
    rhs  = 2*u - (1/2)*uold;
    uold = u;
    u    = H_bdf\rhs;
  end

  % Exact solution: Fourier series for u0 = pi/4
  uex = 0*x;
  for m=1:2:199
    bm = (pi/4) * 2*(1 - cos(m*pi))/(m*pi);
    uex = uex + bm * exp(-nu*(m*pi)^2*time) * sin(m*pi*x);
  end
  err = uex-u;
  emx = max(abs(err))/max(abs(uex));
  el2 = sqrt( (err'*err) / (uex'*uex) );
  elast = enew;

end;

enew  = el2; ratio = elast/enew;

format longe
disp([name '  ' num2str([N nstep el2 ratio ],' %6d  %6d  %e  %e')])
