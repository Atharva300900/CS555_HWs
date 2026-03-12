% q3c_driver.m
% Q3c driver: show improved time convergence
% Plot relative error vs dt for ADI-CN (2nd order) and ADI-BDF3 (3rd order)

close all
clear all
clc

hdr

T = 1.2;

% Pick N large enough that spatial error is small compared to time error
N  = 200;

dts = [0.16000 0.08000 0.04000 0.02000 0.01000 0.00500 0.00250 0.00125];

e_adi  = zeros(size(dts));
e_bdf3 = zeros(size(dts));

fprintf('\n===== Q3c: time convergence at fixed N=%d (T=%.3f) =====\n', N, T);
fprintf('    dt   |  err_ADI-CN   ratio |  err_ADI-BDF3  ratio\n');
fprintf('  ------ | -----------  ----- | ------------  -----\n');

r_adi  = nan(size(dts));
r_bdf3 = nan(size(dts));

for i = 1:length(dts)
    dt = dts(i);

    e_adi(i)  = heat2d_adi_cn(N, dt, T);
    e_bdf3(i) = heat2d_adi_bdf3(N, dt, T);

    if i > 1
        r_adi(i)  = e_adi(i-1)  / e_adi(i);
        r_bdf3(i) = e_bdf3(i-1) / e_bdf3(i);
    end

    fprintf(' %7.5f | %11.2e %6.2f | %12.2e %6.2f\n', ...
            dt, e_adi(i), r_adi(i), e_bdf3(i), r_bdf3(i));
end

% -----------------------------
% Convergence plot: error vs dt
% -----------------------------
figure
loglog(dts, e_adi,  'o-', 'Color','r', 'MarkerSize',10, lw, 2); hold on
loglog(dts, e_bdf3, 's-', 'Color','b', 'MarkerSize',10, lw, 2);
grid on

xlabel('\Delta t', fs, 12)
ylabel('Relative L_2 error', fs, 12)
title('Convergence: ADI-CN (2nd order) vs ADI-BDF3 (3rd order)', fs, 12)

% Reference slope lines O(dt^2) and O(dt^3), anchored near middle dt
j = ceil(length(dts)/2);
dt0 = dts(j);

c2 = e_adi(j)  / (dt0^2);
c3 = e_bdf3(j) / (dt0^3);

ref2 = c2 * (dts.^2);
ref3 = c3 * (dts.^3);

figure('Position', [100 100 600 600]);
loglog(dts, ref2, '--', 'Color','r', lw, 1.5);
loglog(dts, ref3, '--', 'Color','b', lw, 1.5);

legend('ADI-CN', 'ADI-BDF3', 'O(\Delta t^2) ref', 'O(\Delta t^3) ref', ...
       'Location', 'SouthEast')
