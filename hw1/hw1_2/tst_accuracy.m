
% Accuracy comparison: Original CN vs Permuted CN vs ADI
% Fixed dt, varying N (spatial convergence)

dt = 0.01;

fprintf('\n===== Original CN =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1;
N= 20; heat2d
N= 40; heat2d
N= 60; heat2d
N= 80; heat2d
N=100; heat2d
N=120; heat2d
N=160; heat2d
N=200; heat2d

fprintf('\n===== Permuted CN =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1;
N= 20; heat2d_perm
N= 40; heat2d_perm
N= 60; heat2d_perm
N= 80; heat2d_perm
N=100; heat2d_perm
N=120; heat2d_perm
N=160; heat2d_perm
N=200; heat2d_perm
fprintf('\n===== ADI =====\n');
fprintf('   N       dt      nstep    T       error      ratio\n');
e2=1;
N= 20; heat2d_adi
N= 40; heat2d_adi
N= 60; heat2d_adi
N= 80; heat2d_adi
N=100; heat2d_adi
N=120; heat2d_adi
N=160; heat2d_adi
N=200; heat2d_adi
