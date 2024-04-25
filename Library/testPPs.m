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

DrawTrajLVLH3D(state(:, 1:3)*DU, color, linestyle);


end