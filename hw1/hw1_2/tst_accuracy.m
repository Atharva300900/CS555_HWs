
% Accuracy comparison: Original CN vs Permuted CN vs ADI

fprintf('\n===== Original CN =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1; N=80;
dt=.16000; heat2d
dt=.08000; heat2d
dt=.04000; heat2d
dt=.02000; heat2d
dt=.01000; heat2d
dt=.00500; heat2d
dt=.00250; heat2d
dt=.00125; heat2d

fprintf('\n===== Permuted CN =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1; N=80;
dt=.16000; heat2d_perm
dt=.08000; heat2d_perm
dt=.04000; heat2d_perm
dt=.02000; heat2d_perm
dt=.01000; heat2d_perm
dt=.00500; heat2d_perm
dt=.00250; heat2d_perm
dt=.00125; heat2d_perm

fprintf('\n===== ADI =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1; N=80;
dt=.16000; heat2d_adi
dt=.08000; heat2d_adi
dt=.04000; heat2d_adi
dt=.02000; heat2d_adi
dt=.01000; heat2d_adi
dt=.00500; heat2d_adi
dt=.00250; heat2d_adi
dt=.00125; heat2d_adi
