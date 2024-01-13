function COE = MEE2COE(MEE)
% Description: this function converts from Modified Equinoctial Elements to
% Classical Orbital Elements.

if size(MEE, 2) ~= 6
    error('Wrong Dimension of MEE! MEE is expected to be a Nx6 Matrix.')
end

N = size(MEE, 1);
COE = zeros(N, 6);

for i = 1 : N
    % Retrieve MEE from Input
    p = MEE(i, 1);
    l = MEE(i, 2);
    m = MEE(i, 3);
    n = MEE(i, 4);
    s = MEE(i, 5);
    q = MEE(i, 6);
    
    % Compute COE
    e = sqrt(l^2 + m^2);
    a = p/(1 - (l^2 + m^2));
    incl = 2 * atan(sqrt(n^2 + s^2));
    Omega = 2 * atan(s/(n + sqrt(n^2 + s^2)));

    omega = 2 * atan(m/(sqrt(l^2 + m^2) + l)) - Omega;
    omega = fixdomain(omega, -pi, pi);

    nu = q - 2 * atan(m/(l + sqrt(l^2 + m^2)));
    nu = fixdomain(nu, -pi, pi);

    COE(i, :) = [a, e, incl, Omega, omega, nu];

    % % Optional - Check the COE values
    % checkCOE(COE(i, :));

end