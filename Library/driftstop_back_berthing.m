function [value, isterminal, direction] = driftstop_back_berthing(~, Y, exclusion_radius)
% Description: this function identifies the condition when integration must 
% stop, i.e. when natural drift from desired final state has reached a 
% relative distance equal to the exclusion radius.

global DU

if nargin < 3
    exclusion_radius = 15e-3;      % 15 meters
end

% Retrieve Data from State
rho_LVLH = Y(7:9);
rho = norm(rho_LVLH);

value = rho - exclusion_radius/DU;
isterminal = 1;
direction = 0;

end