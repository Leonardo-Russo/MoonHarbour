function MEE = COE2MEE(COE)
% Description: this function computes the Modified Equinoctial Elements
% starting from the Classical Orbital Elements.

if size(COE, 2) ~= 6
    error(['Wrong Dimension of COE!' ...
        'COE is expected to be a Nx6 Matrix.'])
end

N = size(COE, 1);
MEE = zeros(N, 6);

for i = 1 : N

    % Retrieve values from COE
    a = COE(i, 1);
    e = COE(i, 2);
    incl = COE(i, 3);
    Omega = COE(i, 4);
    omega = COE(i, 5);
    nu = COE(i, 6);
    
    % Compute MME
    p = a*(1-e^2);
    l = e*cos(Omega+omega);
    m = e*sin(Omega+omega);
    n = tan(incl/2)*cos(Omega);
    s = tan(incl/2)*sin(Omega);
    q = Omega + omega + nu;

    MEE(i, :) = [p, l, m, n, s, q];

end

end