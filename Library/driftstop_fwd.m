function [value, isterminal, direction] = driftstop_fwd(~, state)
% Description: this function identifies the condition when integration must
% stop, i.e. when natural drift from final controlled state has reached a 
% relative distance of 5m.
% 
% Credits: This function was developed building upon a core reference 
% provided by Dario Sanna.

global DU

rho_vect = state(7:9);
rho = norm(rho_vect);

tol = 1e-6;

value = rho - (5e-3 + tol)/DU;       % condition that stops integration when rho=5m (adimensionalised)
isterminal = 1;
direction = 0;


end