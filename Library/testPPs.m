function testPPs(PPs, tspan)
% Description: this function plots the evolution of the trajectory provided
% in terms of PPs over the time in tspan.

global DU

state = zeros(length(tspan), 3);

for i = 1 : length(tspan)
    state(i, :) = ppsval(PPs, tspan(i));
end

DrawTrajLVLH3D(state(:, 1:3)*DU, '#6efad2', '-.');


end