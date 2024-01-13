function COE = rvPCI2COE(X, mu)
% Description:
% Evaluates the classical orbital parameters @ a given epoch starting 
% from the position and velocity vectors.

% Input: 
% X contains the Planet Centered Inertial position and velocities in a Nx6 shape

% Output:
% COE vector contains the Classical Orbital Elements in the following order
% COE = [a, e, incl, Omega, omega, nu]

N = size(X, 1);
COE = zeros(N, 6);

for i = 1 : N
    
    x = X(i, 1);
    y = X(i, 2);
    z = X(i, 3);
    vx = X(i, 4);
    vy = X(i, 5);
    vz = X(i, 6);
    
    r = sqrt(x^2 + y^2 + z^2);
    rvect = [x, y, z];
    rver = rvect / r;
    
    v = sqrt(vx^2 + vy^2 + vz^2);
    vvect = [vx, vy, vz];
    vver = vvect / v;
    
    hvect = cross(rvect, vvect);
    h = norm(hvect);
    hver = hvect / h;
    
    a = mu / (2 * mu/r - v^2);
    e = sqrt(1 - h^2 / (mu * a));
    
    p = a * (1 - e^2);
    thetaver = cross(hver, rver);
    vr = 1 / r * dot(vvect, rvect);
    
    incl = acos(hver(3));
    Omega = atan2(hver(1) / sin(incl),- hver(2) / sin(incl));
    thetastar = atan2(vr / e * sqrt(p / mu), 1 / e * (p / r - 1));
    thetat = atan2(rver(3) / sin(incl), thetaver(3) / sin(incl));
    omega = wrapToPi(thetat - thetastar);
    
    % % Optional - Check the COE values
    % checkCOE([a, e, incl, Omega, omega, thetastar]);  % Call to your function to check the computed COE values

    COE(i, :) = [a e incl Omega omega thetastar];

end

end
