function showCOE(COEs, idx_ref)
% Description: this function writes in the command window all the classical
% orbital elements at each step.

if nargin < 2
    idx_ref = 0;
end

N = size(COEs, 1);

for i = 1 : N
    fprintf('<strong>Classical Orbital Elements @ %.0f:</strong>\n', i+idx_ref)
    fprintf('a = %.6f\ne = %.6f\ni = %.6f\nOmega = %.6f\nomega = %.6f\nnu = %.6f\n\n', [COEs(i, 1:2), rad2deg(COEs(i, 3:6))])
end

end