function [value, stop] = is_terminal_distance(t, state)
% Assuming y contains position information in the first three elements
% Adjust the indices according to your state vector structure

global DU

rho_vect = state(7:9);
rho = norm(rho_vect);

% Event condition: distance less than 50 meters
value = rho - 50e-3/DU; % Becomes negative when distance is less than 50 m

% Stop integration when distance is less than 50 meters
stop = value <= 0;

end
