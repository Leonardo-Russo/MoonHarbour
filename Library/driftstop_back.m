function [value, isterminal, direction] = driftstop_back(~, Y)
% Description: this function identifies the condition when integration must 
% stop, i.e. when natural drift from desired final state has reached a 
% relative distance of 10m.
% 
% Credits: This function was developed building upon a core reference 
% provided by Dario Sanna.

global DU

% Retrieve Data from State
rho_LVLH = Y(7:9);
rho = norm(rho_LVLH);

% Condition that stops integration when rho = 10m
tol = 10e-3;    % m

value = rho - tol/DU;
isterminal = 1;
direction = 0;

end