function mis = define_misalignment_error(choice, beta, gamma)
% Description: this function describes the misalignment error inside a
% struct.

if nargin < 2
    beta = 0;
end
if nargin < 3
    gamma = 0;
end


mis = struct("name", "Misalignment Parameters");


if choice == "const"
    mis.beta = deg2rad(45);
    mis.gamma = deg2rad(1);
elseif choice == "null"
    mis.beta = 0;
    mis.gamma = 0;
elseif choice == "montecarlo"
    mis.beta = beta;
    mis.gamma = gamma;
else
    error('Unknown Misalignment Error Parameters');
end


end