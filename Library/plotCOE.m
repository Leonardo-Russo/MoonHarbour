function plotCOE(COE, tspan)
% Description: this function creates a plot with the evolution of each COE.

N = length(tspan);

subplot(3, 2, 1)
plot(tspan, COE(:, 1))
xlabel('t')
ylabel('a')

subplot(3, 2, 2)
plot(tspan, COE(:, 2))
xlabel('t')
ylabel('e')

subplot(3, 2, 3)
plot(tspan, rad2deg(COE(:, 3)))
hold on
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, 0*ones(N, 1), 'r--')
xlabel('t')
ylabel('incl')

subplot(3, 2, 4)
plot(tspan, rad2deg(COE(:, 4)))
hold on
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('t')
ylabel('Omega')

subplot(3, 2, 5)
plot(tspan, rad2deg(COE(:, 5)))
hold on
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('t')
ylabel('omega')

subplot(3, 2, 6)
plot(tspan, rad2deg(COE(:, 6)))
hold on
plot(tspan, 180*ones(N, 1), 'r--')
plot(tspan, -180*ones(N, 1), 'r--')
xlabel('t')
ylabel('nu')

end