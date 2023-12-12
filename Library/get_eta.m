function eta = get_eta(MEE)
% Description: this function defines the eta variable useful in the
% propagation of Modified Equinoctial Elements.
% 
% Inputs:
% MEE = [x1, x2, x3, x4, x5, x6]
% 
% Outputs:
% eta = gravitational model help variable "eta"

% Retrieve Data from MEE
x2 = MEE(2);
x3 = MEE(3);
x6 = MEE(6);

eta = 1 + x2 * cos(x6) + x3 * sin(x6);

end