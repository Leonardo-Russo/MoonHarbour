function [value, isterminal, direction] = driftstop_fwd_docking(~, state, gateway_radius)
% Description: this function identifies the condition when integration must
% stop, i.e. when natural drift from final controlled state has reached a 
% relative distance equal to the final distance for docking.

global DU

if nargin < 3
    gateway_radius = 5e-3;
end

rho_vect = state(7:9);
rho = norm(rho_vect);

tol = 1e-6;

value = rho - (gateway_radius + tol)/DU;       % condition that stops integration when rho=5m (adimensionalised)
isterminal = 1;
direction = 0;


end