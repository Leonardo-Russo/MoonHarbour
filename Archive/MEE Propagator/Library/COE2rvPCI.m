function X = COE2rvPCI(COE, mu)
% Description:
% Transformation from Classical Orbital Elements into Planet Centered
% Inertial State.

% Input: 
% COE vector contains the Classical Orbital Elements in the following order
% COE = [a, e, incl, Omega, omega, nu]

% Output:
% X contains the PCI position and velocities in a Nx6 shape

N = size(COE, 1);
X = zeros(N, 6);

for i = 1 : N
    
    % Read values from input
    a = COE(i, 1);
    e = COE(i, 2);
    incl = COE(i, 3);
    Omega = COE(i, 4);
    omega = COE(i, 5);
    nu = COE(i, 6);
    
    r = a * (1 - e^2) / (1 + e * cos(nu));
    x = r * (cos(omega + nu) * cos(Omega) -  cos(incl) * sin(omega + nu) * sin(Omega));
    y = r * (cos(omega + nu) * sin(Omega) +  cos(incl) * sin(omega + nu) * cos(Omega));
    z = r * sin(incl) * sin(omega + nu);
    vx = sqrt(mu / (a * (1 - e^2))) * (- cos(Omega) * (sin(omega + nu) + ...
        e * sin(omega)) - sin(Omega) * cos(incl) * (cos(omega + nu) + e * cos(omega)));
    vy = sqrt(mu / (a * (1 - e^2))) * (cos(Omega) * cos(incl) * (cos(omega + nu) + ...
        e * cos(omega)) - sin(Omega) * (sin(omega + nu) + e * sin(omega)));
    vz = sqrt(mu / (a * (1 - e^2))) * sin(incl) * (cos(omega + nu) + e * cos(omega));

    X(i, :) = [x y z vx vy vz];

end

end
