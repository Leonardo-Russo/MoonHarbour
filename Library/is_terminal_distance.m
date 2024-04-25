function [value, stop] = is_terminal_distance(~, state, terminal_dist)
% Assuming y contains position information in the first three elements
% Adjust the indices according to your state vector structure.

global DU

if nargin < 3
    terminal_dist = 100e-3/DU;
end

rho_vect = state(7:9);
rho = norm(rho_vect);

% Event Condition
value = rho - terminal_dist;    % negative when distance < terminal_dist

% Stop integration when distance is less than terminal_dist
stop = value <= 0;

end
