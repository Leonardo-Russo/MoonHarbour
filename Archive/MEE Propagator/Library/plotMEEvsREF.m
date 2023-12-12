function plotMEEvsREF(MEE, MEE_REF, tspan)
% Description: this function creates a plot with the evolution of each COE.

subplot(3, 2, 1)
plot(tspan, MEE(:, 1))
hold on
plot(tspan, MEE_REF(:, 1))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$p$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 2)
plot(tspan, MEE(:, 2))
hold on
plot(tspan, MEE_REF(:, 2))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$l$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 3)
plot(tspan, MEE(:, 3))
hold on
plot(tspan, MEE_REF(:, 3))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$m$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 4)
plot(tspan, MEE(:, 4))
hold on
plot(tspan, MEE_REF(:, 4))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$n$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 5)
plot(tspan, MEE(:, 5))
hold on
plot(tspan, MEE_REF(:, 5))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$s$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

subplot(3, 2, 6)
plot(tspan, MEE(:, 6))
hold on
plot(tspan, MEE_REF(:, 6))
xlabel('$t \ [days]$', 'Interpreter','latex', 'FontSize', 11)
ylabel('$q$', 'Interpreter','latex', 'FontSize', 11)
legend('Predicted', 'Reference', 'location', 'best')

end