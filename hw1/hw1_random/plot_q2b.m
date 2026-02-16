% Question 2b: Growth factors for various time steppers
% Plot G(z) for z = lambda*dt in [-10, 0]

z = linspace(-10, 0, 500);

G_exact = exp(z);                       % Analytical
G_ef    = 1 + z;                        % Euler Forward
G_eb    = 1 ./ (1 - z);                 % Euler Backward
G_cn    = (1 + z/2) ./ (1 - z/2);      % Crank-Nicolson

figure(1); hold off;
plot(z, G_exact, 'k-',  'linewidth', 2); hold on;
plot(z, G_ef,    'b--', 'linewidth', 2);
plot(z, G_eb,    'r-.', 'linewidth', 2);
plot(z, G_cn,    'm:',  'linewidth', 2);
plot(z, 0*z,     'k-',  'linewidth', 0.5);           % zero line
plot(z, -1+0*z,  'k--', 'linewidth', 0.5);           % G = -1 asymptote

axis([-10 0 -2 1.2]);
xlabel('\lambda \Delta t','FontSize',14);
ylabel('G(\lambda \Delta t)','FontSize',14);
title('Growth Factors','FontSize',15);
legend('Analytical: e^{\lambda\Delta t}', ...
       'Euler Forward: 1 + \lambda\Delta t', ...
       'Euler Backward: 1/(1 - \lambda\Delta t)', ...
       'Crank-Nicolson: (1+z/2)/(1-z/2)', ...
       'Location','SouthEast');
legend boxoff;
grid on;

print('-dpng','-r150','plot_q2b.png');
