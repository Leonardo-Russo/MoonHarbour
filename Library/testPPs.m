function testPPs(PPs, tspan, color, linestyle)
% Description: this function plots the evolution of the trajectory provided
% in terms of PPs over the time in tspan.

if nargin < 3
    color = '#14dde0';
end
if nargin < 4
    linestyle = '-.';
end

global DU

state = zeros(length(tspan), 3);

for i = 1 : length(tspan)
    state(i, :) = ppsval(PPs, tspan(i));
end

x_max = max(state(:, 1))*DU;
x_min = min(state(:, 1))*DU;
y_max = max(state(:, 2))*DU;
y_min = min(state(:, 2))*DU;

fig = figure('name', 'testpps');
DrawTrajLVLH3D(state(:, 1:3)*DU, color, linestyle);
view(0, 90)
xlim([x_min, x_max]);
ylim([y_min, y_max]);
res = '-r600';
print(fig, 'Output/reference traj.png', '-dpng', res);


end