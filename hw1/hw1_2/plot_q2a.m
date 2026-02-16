% Plot for 2a: Sparsity patterns of original vs permuted Cholesky factors
% Matrix size: 2500 x 2500  (N=51, m=50, m^2=2500)

hdr
N = 51; m = N-1;
Lx=1; Ly=1; dt=0.01; nu=0.05;

dx=Lx/N; e = ones(m,1); A = spdiags([-e 2*e -e],-1:1,m,m);
Ax = nu*A/(dx*dx); Ix = speye(m);
dy=Ly/N; A = spdiags([-e 2*e -e],-1:1,m,m);
Ay = nu*A/(dy*dy); Iy = speye(m);

HL = kron(Iy,Ix) + (dt/2)*(kron(Iy,Ax)+kron(Ay,Ix));

fprintf('Matrix size: %d x %d\n', size(HL));

% Original Cholesky
L_orig = chol(HL,'lower');

% Permuted Cholesky
p = symamd(HL);
L_perm = chol(HL(p,p),'lower');

fprintf('Original nnz(L): %d\n', nnz(L_orig));
fprintf('Permuted nnz(L): %d\n', nnz(L_perm));
fprintf('Ratio: %.2fx sparser\n', nnz(L_orig)/nnz(L_perm));

figure(1);
subplot(1,2,1);
spy(L_orig);
title(sprintf('Original Cholesky L\nnnz = %d', nnz(L_orig)),'fontsize',13);

subplot(1,2,2);
spy(L_perm);
title(sprintf('Permuted Cholesky L\nnnz = %d', nnz(L_perm)),'fontsize',13);

print('-dpng','-r150','plot_q2a.png');
