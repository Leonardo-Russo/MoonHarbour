function plotCOEvsREF(COE, COE_REF, tspan)
% Description: this function creates a plot with the evolution of each COE.

N = length(tspan);

subplot(3, 2, 1)
plot(tspan, COE(:, 1))
hold on
plot(tspan, COE_REF(:, 1))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$a$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 2)
plot(tspan, COE(:, 2))
hold on
plot(tspan, COE_REF(:, 2))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$e$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 3)
plot(tspan, rad2deg(COE(:, 3)))
hold on
plot(tspan, rad2deg(COE_REF(:, 3)))
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, 0*ones(N, 1), 'r--')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$i$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 4)
plot(tspan, rad2deg(COE(:, 4)))
hold on
plot(tspan, rad2deg(COE_REF(:, 4)))
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$\Omega$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 5)
plot(tspan, rad2deg(COE(:, 5)))
hold on
plot(tspan, rad2deg(COE_REF(:, 5)))
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$\omega$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 6)
plot(tspan, rad2deg(COE(:, 6)))
hold on
plot(tspan, rad2deg(COE_REF(:, 6)))
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$\theta^*$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')


end